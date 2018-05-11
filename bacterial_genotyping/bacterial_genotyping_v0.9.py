#!/usr/bin/env python

'''
By Sean Lovett
8/18/2016

This pipeline performs bacterial genotyping. It automates the process of screening for repetitive and low complexity sequence
using Prinseq and custom scripts, trimming adapters with Cutadapt, quality trimming/filtering with Prinseq, alignment with
Bowtie2, alignment filtering/sorting with Samtools, duplicate removal and adding of readgroup information with Picard,
calculation of depth of coverage with Samtools, plotting of depth of coverage with custom scripts, merging of bam files with
Samtools, variant calling with GATK's UnifiedGenotyper, variant annotation with SnpEff, and filtering and summarization of
annotated variants with custom scripts.

The way to call most programs can be changed from the commandline.
Scripts written at RIID are expected to be in the directory corresponding to the --script_dir flag.

Not tested yet:
  Single sample mode.
  
Future features:
  More options - especially for coverage plots.
  gzipped input
  Incorporate logic from external scripts instead of calling them from the commandline

Necessary inputs:
  Reference fasta (edit headers to only include accession number)
  GFF file (accession must match ref fasta)
  Reads (must be paired-end)
  Tab-delimited text file of format (sample\tR1.fastq\tR2.fastq)
Outputs:
  Merged bam with separate read group for each sample
  (Optional) Depth report for each sample and coverage plots for each sample
  Raw vcf file from UnifiedGenotyper
  Annotated and summarized vcf file (suffix = .ann.vcf.summary.final)
'''

from __future__ import print_function
from subprocess import Popen, PIPE
import optparse, time, sys

def main():
  #To parse command line
  usage = "usage: %prog -i input.txt -r ref.fasta -g ref.gff [options] > commandlog.txt"
  p = optparse.OptionParser(usage)
  
  input_opts = optparse.OptionGroup(p, 'Input files')
  input_opts.add_option('-i', '--input', help='Input read information. Tab-delimited, 3 columns per sample, samplename R1.fastq R2.fastq [None,REQD]')
  input_opts.add_option('-r', '--ref', help='Input reference fasta (edit headers to only include accession number [None,REQD]')
  input_opts.add_option('-g', '--gff', help='Input reference GFF. Must use even if snpEff db is predefined. Accession(s) must match ref fasta [None, REQD]')
  p.add_option_group(input_opts)

  func_opts = optparse.OptionGroup(p,'Optional functions')
  func_opts.add_option('--depth', default=False, action="store_true", help='Generate depth report using samtools [False]')
  func_opts.add_option('--clean', default=False, action="store_true", help='Clean up intermediate files, leaving only merged bam, final vcf, and depth reports [False]')
  func_opts.add_option('--dry_run', default=False, action="store_true", help="Dry run. Only print commands [False]")
  p.add_option_group(func_opts)
  
  param_opts = optparse.OptionGroup(p, 'Parameters')
  param_opts.add_option('-t', '--threads', default=2, type="int", help="Number of threads to use [2]")
  param_opts.add_option('-w', '--window', default=50, type="int", help="Window size for cut_ref_fastq.py [50]")
  param_opts.add_option('--scorediffmax', type="int", default=0, help="Score difference between AS and XS must be be less than this for a region not to be repetitive [0]")
  param_opts.add_option('--lcthresh', type="int", default=3, help="Low complexity threshold for marking with prinseq [3]")
  param_opts.add_option('-m', '--minlen', default=100, type="int", help="Minimum read length after filtering [100]")
  param_opts.add_option('--mapq', type="int", default=10, help="Minimum mapping quality to keep read [10]")
  param_opts.add_option('--use_db', help="snpEff db to use. One will be created from gff if not set [None]")
  param_opts.add_option('--min_freq', type="int", default=10,help="Minimum percent frequency to report variant [10]")
  param_opts.add_option('--min_depth', type="int", default=20, help="Minimum depth to report variant [20]")
  p.add_option_group(param_opts)
  
  prog_opts = optparse.OptionGroup(p, 'Program paths')
  prog_opts.add_option('--java', default = "/usr/lib/jvm/java-1.7.0/jre/bin/java", help="How to call java. Must be 1.7+ [/usr/lib/jvm/java-1.7.0/jre/bin/java]")
  prog_opts.add_option('--cut', default = "/data/_scripts/slovett/cutadapt_paired_nosispa_2.0.py", help="How to call cutadapt [/data/_scripts/slovett/cutadapt_paired_nosispa_2.0.py]")
  prog_opts.add_option('--prinseq', default="/data/_pipeline_software/prinseq-lite-0.20.3/prinseq-lite.pl", help='How to call prinseq [/data/_pipeline_software/prinseq-lite-0.20.3/prinseq-lite.pl]')
  prog_opts.add_option('--bowtie', default="/data/_pipeline_software/bowtie2-2.1.0/bowtie2", help="How to call Bowtie2 [/data/_pipeline_software/bowtie2-2.1.0/bowtie2]")
  prog_opts.add_option('--picard', default="/data/_pipeline_software/picard-tools/dist/picard.jar", help='How to call picard [/data/_pipeline_software/picard-tools/dist/picard.jar]')
  prog_opts.add_option('--gatk', default="/data/_pipeline_software/GenomeAnalysisTK.jar", help="How to call GATK [/data/_pipeline_software/GenomeAnalysisTK.jar]")
  prog_opts.add_option('--snpeff', default="/data/_pipeline_software/snpEff", help="snpEff basedir. Directory in which snpEff.jar resides, NOT THE FILE ITSELF. [/data/_pipeline_software/snpEff]")
  prog_opts.add_option('--samtools', default="/data/_pipeline_software/samtools-1.3/samtools", help="How to call samtools, must be 1.3+ [/data/_pipeline_software/samtools-1.3/samtools]")
  prog_opts.add_option('--script_dir', default="/data/_scripts/slovett/bacterial_genotyping", help="Directory containing needed scripts [/data/_scripts/slovett]")
  p.add_option_group(prog_opts)
  
  out_opts = optparse.OptionGroup(p, 'Output files')
  out_opts.add_option('-b', '--bam', default="merged.bam", help="Name for output merged bam [merged.bam]")
  p.add_option_group(out_opts)

  opts, args = p.parse_args()

  print("# COMMAND: %s\n" % " ".join(sys.argv))

  # Handle obviously bad command lines
  if opts.input is None: print("Missing input read info.")
  if opts.ref is None: print("Missing input reference fasta.")
  if opts.gff is None: print("Missing input reference GFF")
  if opts.input is None or opts.ref is None or opts.gff is None: exit(-1)
  
  reads = {}
  bams = []
  depths = []
  files_to_remove = []

  # Read input file and create dict of tuples key:value = sample_name:(read1,read2)
  with open(opts.input) as in_reads:
    for line in in_reads:
      cols = line.strip().split("\t")
      reads[cols[0]] = (cols[1], cols[2])

  # Prep reference fasta
  print("# Preparing reference files")
  run_faidx(opts)
  run_dict(opts)
  ref_base = run_btindex(opts)
  print("")

  # Get info about repetetive and low complexity regions
  print("# Scanning for repetitive and low-complexity sequence in the reference")
  cut_ref_fastq = run_cut_ref(opts)
  rep_regions_bed = run_rep_regions(opts, cut_ref_fastq)
  lc_regions_bed = run_lc_regions(opts, cut_ref_fastq)
  files_to_remove.extend((cut_ref_fastq, lc_regions_bed))
  print("")

  # Process each sample individually, generating one bam per sample
  for sample in sorted(reads):
    print("# Processing %s" % sample)
    read1_cut, read2_cut = run_cutadapt(sample, reads[sample][0], reads[sample][1], opts)
    read1_ps, read2_ps = run_prinseq(sample, read1_cut, read2_cut, opts)
   
    bam = run_align(sample, ref_base, read1_ps, read2_ps, opts)
    rmdups = run_rmdups(sample, bam, opts)
    bam_rg = run_addrg(sample, rmdups, opts)

    if opts.depth:
      depth = run_depth(bam_rg, opts)
      depths.append(depth)
   
    bams.append(bam_rg)
    print("")
    # files_to_remove.extend((read1_cut, read2_cut, read1_ps, read2_ps, bam, rmdups, bam_rg))
    files_to_remove.extend((read1_cut, read2_cut, bam, rmdups))

  if opts.depth:
    print("# Plotting depths of coverage")
    if opts.dry_run:
      print("# These commands are not printed during a dry run")
    depth_plots = run_depth_plots(depths, opts)
    print("")

  # Merge bams and call variants
  print("# Merging bams and calling variants")
  if len(bams) > 1:
    single = False
    run_merge(bams, opts)
    vcf_file = call_vars(opts.bam, opts)
  else:
    single = True
    run_indexbam(bams[0])
    vcf_file = call_vars(bams[0],opts)
  print("")

  # Run SNPEff
  # Would be nice to check for existence of db and not attempt to create one if it already exists
  print("# Running SNPEff")
  if not opts.use_db: opts.use_db = make_snpeff_db(opts)
  snpeff_vcf = run_snpeff(opts.use_db, vcf_file, opts)
  print("")

  # Summarize and filter snpeff output
  print("# Generating final output")
  summary_vcf = run_summary(snpeff_vcf, single, opts)

  # Generate lc and rep annotations  
  lc_annot = run_bed_to_annot(lc_regions_bed, summary_vcf, opts, "lc")
  rep_annot = run_bed_to_annot(rep_regions_bed, summary_vcf, opts, "rep")

  files_to_remove.extend((lc_annot, rep_annot))

  # Merge vcf with lc and rep annotations to create final output  
  if not opts.dry_run: final_vcf = add_lc_rep_annots(summary_vcf, lc_annot, rep_annot)
  print("")

  if opts.clean:
    print("# Cleaning up")
    cleanup(files_to_remove, opts)
    print("")
    
#--------------------------------------------------------------------------------------------------

def run_depth_plots(depth_files, opts):
  chroms = set()

  if not opts.dry_run:
    with open(depth_files[0], "r") as fin:
      for line in fin:
        cols = line.strip().split("\t")
        chroms.add(cols[0])

  depth_dict = {}
  for chrom in chroms:
    depth_dict[chrom] = []
    for depth in depth_files:
      out_file = "%s_%s.depth" % (".".join(depth.split(".")[:-1]),chrom)
      depth_dict[chrom].append(out_file)
      grep_cline= "grep %s %s > %s" % (chrom, depth, out_file)
      print(grep_cline)
      
      if not opts.dry_run:
        run_grep = Popen(grep_cline, shell=True, stdout=PIPE, stderr=PIPE)
        run_grep.wait()

  cov_fofns = []        
  for chrom, depth_list in depth_dict.iteritems():
    cov_fofn = "%s.depths" % chrom
    cov_fofns.append(cov_fofn)
    with open(cov_fofn, "w") as fout:
      for depth in depth_list:
        fout.write("%s\n" % depth)
        
  for fofn in cov_fofns:
    covplot_cline = ("Rscript %s/covplots_cli.R %s 5000 2 %s" % (opts.script_dir, fofn,".".join(fofn.split(".")[:-1])))
    print(covplot_cline)
    
    if not opts.dry_run:
      run_covplot = Popen(covplot_cline, shell=True, stdout=PIPE, stderr=PIPE)
      run_covplot.wait()
  
  
def cleanup(files_to_remove, opts):
  cline = "rm %s *.bai *.dict *.fai *singletons.fastq *bt2 *.sam snpEff_genes.txt snpEff_summary.html" % " ".join(files_to_remove)
  print(cline) 

  if not opts.dry_run:
    Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)


def add_lc_rep_annots(vcf, lc_annot, rep_annot):
  out_name = "%s.final" % vcf

  with open(vcf,"r") as vcf_in:
    with open(out_name,"w") as vcf_out:
      i=0
      with open(lc_annot, "r") as lc_in: lc_annots = lc_in.readlines()
      with open(rep_annot, "r") as rep_in: rep_annots = rep_in.readlines()
      for line in vcf_in:
        if line.startswith("#"):
          vcf_out.write(line)
        elif line.startswith("CHROM"):
          cols = line.strip().split("\t")
          vcf_out.write("%s\tLC\tREP\t%s\n" % ("\t".join(cols[0:2]), "\t".join(cols[2:])))
        else:
          vcf_out.write(lc_annots[i].strip())
          rep_cols = rep_annots[i].strip().split("\t")
          vcf_out.write("\t%s\t" % rep_cols[2])
          vcf_cols = line.split("\t")
          vcf_out.write("\t".join(vcf_cols[2:]))
          i+=1
  return out_name


def run_depth(bam, opts):
  out_name = "%s.depth" % ".".join(bam.split(".")[:-1])
  cline = "%s depth -a %s > %s" % (opts.samtools, bam, out_name)
  print(cline)

  if opts.dry_run:
    return out_name
  else:
    run_samdepth=Popen(cline,shell=True,stdout=PIPE,stderr=PIPE)
    run_samdepth.wait()
  
  return out_name


def run_summary(vcf, single, opts):
  out_name = "%s.summary" % vcf
  if single:
    cline="%s/summarize_vcf_v2.py -i %s -m %s -d %s --single > %s" % (opts.script_dir, vcf, opts.min_freq, opts.min_depth, out_name)
  else:
    cline="%s/summarize_vcf_v2.py -i %s -m %s -d %s > %s" % (opts.script_dir, vcf, opts.min_freq, opts.min_depth, out_name)
  print(cline)

  if opts.dry_run:
    return out_name
  else:  
    run_summarize=Popen(cline,shell=True,stdout=PIPE,stderr=PIPE)
    run_summarize.wait()
  
  return out_name


def run_bed_to_annot(bed, vcf, opts, prefix):
  out_name = "%s_annot.txt" % prefix
  cline = "%s/find_rep_snps.py -i %s -l %s -o %s -w %s" % (opts.script_dir, bed, vcf, out_name, opts.window)
  print(cline)

  if opts.dry_run:
    return out_name
  else:  
    run_bedtoannot=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_bedtoannot.wait()
  
  return out_name


def run_lc_regions(opts, cut_ref_fastq):
  cline = "%s/find_lc_regions.py -i %s -p %s -o lc_regions.bed -t %s" % (opts.script_dir, cut_ref_fastq, opts.prinseq, opts.lcthresh)
  print(cline)
  
  if opts.dry_run:
    return "lc_regions.bed"
  else:
    run_lcregions=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_lcregions.wait()
  
  return "lc_regions.bed"


def run_rep_regions(opts, cut_ref_fastq):
  cline = "%s/rep_regions_bybase.py -q %s -r %s -g %s -o rep_regions -p %s -d %s" % (opts.script_dir, cut_ref_fastq, opts.ref, opts.gff, opts.threads, opts.scorediffmax)
  print(cline)

  if opts.dry_run:
    return "rep_regions_d%s_repregions.bed" % opts.scorediffmax
  else:
    run_repregions=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_repregions.wait()

  return "rep_regions_d%s_repregions.bed" % opts.scorediffmax


def run_cut_ref(opts):
  out_name = "%s.fastq" % ".".join(opts.ref.split(".")[:-1])
  cline = "%s/cut_ref_fastq.py -f %s -o %s -w %s" % (opts.script_dir, opts.ref, out_name, opts.window)
  print(cline)

  if opts.dry_run:
    return out_name
  else:  
    run_cutref=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_cutref.wait()

  return out_name


def run_snpeff(db,vcf,opts):
  out_name = "%s.ann.vcf" % ".".join(vcf.split(".")[:-1])
  cline = "%s -Xmx4g -jar %s/snpEff.jar eff -ud 10 -o vcf -c %s/snpEff.config %s %s > %s" % (opts.java, opts.snpeff,opts.snpeff,db,vcf,out_name)
  print(cline)

  if opts.dry_run:
    return out_name
  else:  
    run_snpeff=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_snpeff.wait()

  return out_name


def make_snpeff_db(opts):
  db_name = ".".join(opts.ref.split(".")[:-1])
  
  if db_name in open("%s/snpEff.config" % opts.snpeff).read():
    print("# Using existing database with name %s" % db_name)
    return db_name

  if opts.dry_run:
    return db_name
  else:
    mkdir_cline="mkdir %s/data/%s" % (opts.snpeff,db_name)
    print(mkdir_cline)
    
    run_mkdir=Popen(mkdir_cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_mkdir.wait()
    
    cp_ref_cline="cp %s %s/data/%s/sequences.fa" % (opts.ref, opts.snpeff, db_name)
    print(cp_ref_cline)
    
    run_cp_ref=Popen(cp_ref_cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_cp_ref.wait()
    
    cp_gff_cline="cp %s %s/data/%s/genes.gff" % (opts.gff, opts.snpeff, db_name)
    print(cp_gff_cline)
    
    run_gff_ref=Popen(cp_gff_cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_gff_ref.wait()
    
    with open("%s/snpEff.config" % opts.snpeff, "a") as config_file:
      config_file.write("\n%s.genome : %s\n" % (db_name, db_name))
    
    build_db_cline="%s -jar %s/snpEff.jar build -gff3 -v %s" % (opts.java, opts.snpeff, db_name)
    print(build_db_cline)
    
    run_build_db=Popen(build_db_cline, shell=True, stdout=PIPE, stderr=PIPE)
    run_build_db.wait()

  return db_name


def call_vars(bam, opts):
  out_name = "%s_variants.vcf" % ".".join(bam.split(".")[:-1])
  cline = "%s -Xmx32g -jar %s -T UnifiedGenotyper -R %s -I %s -o %s -ploidy 1 -glm both -nt %s" % (opts.java, opts.gatk, opts.ref, bam, out_name, opts.threads)
  print(cline)

  if opts.dry_run:
    return out_name
  else:  
    gatk_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    gatk_run.wait()
  
  return out_name


def run_dict(opts):
  cline = "%s -jar %s CreateSequenceDictionary R=%s O=%s.dict" % (opts.java, opts.picard, opts.ref, ".".join(opts.ref.split(".")[:-1]))
  print(cline)

  if opts.dry_run:
    return
  else:
    dict_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    dict_run.wait()


def run_faidx(opts):
  cline="%s faidx %s" % (opts.samtools, opts.ref)
  print(cline)

  if opts.dry_run:
    return
  else:
    faidx_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    faidx_run.wait()


def run_indexbam(bam, opts):
  index_cline = "%s index %s" % (opts.samtools, bam)
  print(index_cline)


  if opts.dry_run:
    return
  else:  
    index_run=Popen(index_cline, shell=True, stdout=PIPE, stderr=PIPE)
    index_run.wait()


def run_merge(bams, opts):
  merge_cline = "%s merge -@ %s %s %s" % (opts.samtools, opts.threads, opts.bam, " ".join(bams))
  print(merge_cline)

  if opts.dry_run:
    return
  else:    
    merge_run=Popen(merge_cline, shell=True, stdout=PIPE, stderr=PIPE)
    merge_run.wait()
  
  run_indexbam(opts.bam, opts)

 
def run_addrg(sample,bam,opts):
  out_name = "%s_readgroups.bam" % ".".join(bam.split(".")[:-1])
  cline = "%s -jar %s AddOrReplaceReadGroups I=%s LB=%s ID=%s PL=illumina SM=%s PU=%s O=%s" % (opts.java, opts.picard, bam, sample, sample, sample, time.strftime("%Y%m%d"), out_name)
  print(cline)

  if opts.dry_run:
    return out_name
  else:    
    addrg_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    addrg_run.wait()
  
  return out_name


def run_rmdups(sample, bam, opts):
  out_name = "%s_rmdups.bam" % ".".join(bam.split(".")[:-1])
  metrics = "%s_rmdups.metrics" % ".".join(bam.split(".")[:-1])
  cline = "%s -XX:ParallelGCThreads=%s -jar %s MarkDuplicates REMOVE_DUPLICATES=true ASSUME_SORTED=true O=%s I=%s M=%s" % (opts.java, opts.threads, opts.picard, out_name, bam, metrics)
  print(cline)

  if opts.dry_run:
    return out_name
  else:    
    rmdups_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    rmdups_run.wait()
  
  return out_name


def run_btindex(opts):
  refindex_base = ".".join(opts.ref.split(".")[:-1])
  refindex_cline = "%s-build %s %s" % (opts.bowtie, opts.ref, refindex_base)
  print(refindex_cline)

  if opts.dry_run:
    return refindex_base
  else:    
    refindex_run=Popen(refindex_cline, shell=True, stdout=PIPE, stderr=PIPE)
    refindex_run.wait()
  
  return refindex_base


def run_align(sample, base, read1, read2, opts):
  out_name = "%s.bam" % sample
  bowtie_cline = "%s -x %s -1 %s -2 %s -p %s -X 1000 2> %s_bowtie.out | %s view -q %s -f 3 -bS - | %s sort -@ %s -o %s -" % (opts.bowtie, base, read1, read2, opts.threads, sample, opts.samtools, opts.mapq, opts.samtools, opts.threads, out_name)
  print(bowtie_cline)

  if opts.dry_run:
    return out_name
  else:    
    bowtie_run=Popen(bowtie_cline, shell=True, stdout=PIPE, stderr=PIPE)
    bowtie_run.wait()
  
  return out_name


def run_prinseq(sample, read1, read2, opts):
  out_name = "%s_cutadapt_prinseq" % (sample)
  cline = "%s -fastq %s -fastq2 %s -out_bad null -out_good %s -min_len %s -trim_qual_right 30 -trim_qual_type min -trim_qual_window 5 -min_qual_mean 20 2> %s_prinseq.err" % (opts.prinseq, read1, read2, out_name, opts.minlen, sample)
  print(cline)

  read1_ps = "%s_1.fastq" % out_name
  read2_ps = "%s_2.fastq" % out_name
  
  if opts.dry_run:
    return read1_ps, read2_ps
  else:  
    prinseq_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    prinseq_run.wait()
  

  
  return read1_ps, read2_ps


def run_cutadapt(samplename, read1, read2, opts):
  cline = "%s -1 %s -2 %s -l %s" % (opts.cut, read1, read2, opts.minlen)
  print(cline)

  read1_cut = '%s_cutadapt_paired.fastq' % '.'.join(read1.split('.')[:-1])
  read2_cut = '%s_cutadapt_paired.fastq' % '.'.join(read2.split('.')[:-1])

  if opts.dry_run:
    return read1_cut, read2_cut
  else:
    cut_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
    cut_run.wait()
  
  return read1_cut, read2_cut

  
if __name__ == "__main__":
  main()
