#!/usr/bin/env python

from __future__ import division
import optparse, os
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def main():
    usage = '%prog [options] fastq1 [fastq2 ...]'
    p = optparse.OptionParser()
    #Inputs
    p.add_option('-1', '--one', help='Read 1 fastq file. Can be gzipped or not [None, REQD]')
    p.add_option('-2', '--two', help='Read 2 fastq file. Can be gzipped or not [None, REQD]')
    
    #How to call programs
    p.add_option('--ca', default='cutadapt1.2', help='How to call cutadapt. Need to use version >=1.2 [cutadapt1.2]')
    p.add_option('--pseq', default='prinseq-lite_0.2.2.pl', help='How to call prinseq. [prinseq-lite_0.2.2.pl]')
    p.add_option('--ray', default='Ray2.2', help='How to call Ray. [Ray2.2]')
    p.add_option('--cap', default='cap3', help='How to call cap3. [cap3]')
    p.add_option('--bb', default='bowtie_buildout_latest_wQualchecks.py', help='How to call bowtie_buildout. [bowtie_buildout_latest_wQualchecks.py]')
    p.add_option('--blast', default='split_blasts_multi_faster.py', help='How to call blast function. [split_blasts_multi_faster.py]')
    
    #Cutadapt Settings
    p.add_option('--i1', default='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', help="Illumina adaptor to be trimmed from the 3' end of R1. [GATCGGAAGAGCACACGTCTGAACTCCAGTCAC]")
    p.add_option('--i2', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT', help='Reverse complement of the Illumina adaptor that will be trimmed from the 3; end of R2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')
    p.add_option('-a', '--sispa', default='GCCGGAGCTCTGCAGATATC', help="Sispa adaptor to use. Will be trimmed from 5' end, the reverse complement will be trimmed from the 3' end Ventor=GCCGGAGCTCTGCAGATATC; Columbia=CGCCGTTTCCCAGTAGGTCTC [GCCGGAGCTCTGCAGATATC]")
    
    #Prinseq Settings
    p.add_option('-l', '--minLength', default=50, type='int', help='Minimum length to keep (after primer clipping) [50]')
    p.add_option('-q', '--meanQual', default=20, type='int', help='Minimum mean quality for a read to be kept [20]')
    p.add_option('--derep', default='14', help='derep setting for prinseq. Use the --NOderep flag to turn replicate removal off [14]')
    p.add_option('--NOderep', default=False, action='store_true',  help='Use this flag if you do not want to remove duplicates [False]')
    p.add_option('--lcMethod', default='dust',  help='Tells prinseq which method to use to get rid of low complexity sequences [dust]')
    p.add_option('--lcThresh', default='3',  help='Low complexity threshold [3]')
    p.add_option('--trimRight', default=20, type='int',  help='All bases below this quality will be trimmed from the right of the sequences. Must be an integer [20]')

    #Ray Settings
    p.add_option('-n', '--procs', default=10, type='int', help='Number of processors used for de novo assembly. Same # procs will be used for the blasting & buildout that follows [10]')
    p.add_option('-k', '--kmer', default=21, type='int', help='kmer size to use for de novo assembly [21]')
    p.add_option('--rayOut', help='Output directory for Ray [rayout_k##]')
    
    #Blasting to sequences expected to be similar to the contigs you are building
    p.add_option('--ns', '--nucSubject', help='Fasta file of nucleotide sequences to compare the query sequences to. Will format automatically. [None]')
    p.add_option('--ps', '--protSubject', help='Fasta file of protein sequences to compare the query sequences to. Will format automatically. [None]')
    p.add_option('--temp', default='./temp', help='Name for temporary working directory. This will be created at the begining of the script and then removed at the end. Same name will be used for bowtie buildout [./temp]')
    p.add_option('-b', '--blastType', help='Type of blast to run. Options "blastn", "blastx", "blastp", "tblastx", "tblastn" [blastn or blastx or blastn,blastx]')
    p.add_option('--dontIndex', default=False, action='store_true',  help="Use this flag if you don't want the script to try to index the database. This is necessary for complex databases like nt and nr")
    p.add_option('--task', default='megablast,dc-megablast,blastn', help='Type of blastn to run. Options "blastn", "dc-megablast", "megablast" [megablast,dc-megablast,blastn]')
    p.add_option('--evalue', default='1e-20', help='Maximum evalue for hit to be recorded [1e-20]')
    p.add_option('--numHits', type='int', default=5, help='Integer specifying the number of blast hits to report per query. [5]')
    p.add_option('--numHsps', type='int', default=1, help='Integer specifying the number of alignments to report per query/subject pair. [1]')
    p.add_option('--minContigLen', default=150, type='int', help='Minimum length for a contig to be utilized in the bowtie buildout [150]')
    
    #Bowtie buildout settings
    p.add_option('--bbOut', default='buildout', help='Output directory for bowtie_buildout [buildout]')
    p.add_option('--covThresh', type='int', default=3, help='Minimum level of coverage required to change consensus. [3]')
    p.add_option('--baseQual', type='int', default=20, help='Minimum base quality for a base to be used in the evalution of the consensus. [20]')
    p.add_option('--finalCov', type='int', default=5, help='Minimum level of coverage required to keep a base in the final assembly. If a base is below this threshold, it is changed to an "N" [5]')
    p.add_option('--finalMapQ', type='int', default=20, help='Minimum mapping quality for a read to be used in the pileup generation for the final quality check [20]')
#    p.add_option('--finalMaxCov', type='int', default=500, help='Max per base coverage to be used in the pileup generation for the final quality check [500]')
#    p.add_option('--Nfinal', type='int', default=15, help='# of Ns to add to the beginning and end of each contig prior to the final mapping. This allows for easier mapping at the ends of the contigs.[15]')
#    p.add_option('--propThresh', type='float', default=0.5, help='Requires that a greater proportion of the reads (than this threshold) must support an alternative for the consensus to be changed. [0.5]')
#    p.add_option('--offset', type='int', default=33, help='Base quality offset used in the pileup. I believe the default in samtools is Sanger(33) [33]')
#    p.add_option('--NoReqSameBase', default=False, action='store_true', help='Use this flag to turn off the default behavior, which requires at least --finalCov coverage of a single base to be in the final assembly')
    
    
    #Optional starting places
    #Can currently start in two intermediate places, after cutadapt (bbR1 and bbR2) or after blasting Ray contigs with conInt, in addition to bbR1 and bbR2.
    #Testing other options
    p.add_option('--conInt', help='Use this flag if you already have generated your contigs of interest. If used, you must also provide the reads you want to buiild with using --bbR1 & --bbR2 [None]')
    p.add_option('--bbR1', help='Read 1 fastq file to be used in bowtie buildout. This only neeeds to be provided if you are using an optional starting place')
    p.add_option('--bbR2', help='Read 2 fastq file to be used in bowtie buildout. This only neeeds to be provided if you are using an optional starting place') 
    p.add_option('--pspR1', help='Paired R1 file from prinseq')
    p.add_option('--pspR2', help='Paired R2 file from prinseq')
    p.add_option('--pssR1', help='Singleton R1 file from prinseq')
    p.add_option('--pssR2', help='Singleton R2 file from prinseq')
    p.add_option('--rayDone', default=False, action='store_true',  help="Use this flag if Ray assembly has already been done")

    
    opts, args = p.parse_args()
    
    if not opts.conInt:
        if not opts.rayOut: opts.rayOut = 'rayout_k%d' % opts.kmer
    
        #Set default blast type if none provided
        if not opts.blastType:
            if opts.ns and opts.ps: opts.blastType='blastn,blastx'
            elif opts.ns: opts.blastType='blastn'
            elif opts.ps: opts.blastType='blastx'
            else: '!!!ERROR: must provide at least one subject fasta!!!!'
        
        #Only run this if no paired adaptor clipped files are provided
        if not opts.bbR1 or not opts.bbR2:
            #Store reverse complement of the sispa primer
            opts.sispa_revcomp = rev_comp(opts.sispa)
    
            #Run cutadapt
            one_cut = run_cutadapt(opts.one, opts.i1, opts)
            two_cut = run_cutadapt(opts.two, opts.i2, opts)

            #Make new fastqs with just reads that are still paired. Delete intermediate files
            one_cut_paired, two_cut_paired = paired_sort(one_cut, two_cut)
            delete_files(one_cut, two_cut)	
        #If the user provided paired adaptor clipped files, this will assign these to the correct variables
        else:
            one_cut_paired = opts.bbR1
            two_cut_paired = opts.bbR2
        
        #Prinseq is only run if the output files from prinseq are not provided by the user
        if not opts.pspR1 or not opts.pspR2 or not opts.pssR1 or not opts.pssR2:
            #Run everything through prinseq
            base_good_name = run_prinseq(one_cut_paired, two_cut_paired, opts)
        else:
            base_good_name = '%s_good' % opts.pspR1.split('_good')[0]

        #Ray is not run if the flag --rayDone is not thrown
        if not opts.rayDone:
            #Run Ray assembly
            run_ray2(base_good_name, opts)
    
        #Blast contigs against the reference sequences that were provided to get contigs of interest
        cons_interest = run_blast(opts)
    
    else: 
        cons_interest = opts.conInt
        one_cut_paired = opts.bbR1
        two_cut_paired = opts.bbR2
        
    #Run bowtie buildout to try to extend the contigs using the cutadapt-only fastas
    #Currently I'm running cap3 prior to building out
    if cons_interest:
        cons_interest = os.path.abspath(cons_interest)
        to_build_fasta = run_cap3(cons_interest, opts)
        run_buildout(to_build_fasta, one_cut_paired, two_cut_paired, opts)
    else: print 'No contigs built that match the provided reference'
#-------------------------------End of main()------------------------

def run_buildout(start_fasta, r1, r2, opts):
    cmd='%s -1 %s -2 %s -r %s -o %s --procs %d --temp %s --baseQual %d --covThresh %d --finalMapQ %d --finalCov %d >>bb_stdout &>>bb_stderr' % (opts.bb, r1, r2, start_fasta, opts.bbOut, opts.procs, opts.temp, opts.baseQual, opts.covThresh, opts.finalMapQ, opts.finalCov)
    print cmd
    bbit=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    bbit.wait()

def run_cap3(fasta, opts):
    cmd='%s %s' % (opts.cap, fasta)
    cap3_join=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    cap3_join.wait()
    
    if cap3_improved(fasta):
        #Join CAP3 'conitgs' and 'singlets' into a new fasta file, and save this as the current reference
        new_fasta=join_cap3_contigs(fasta)
        clean_cap3(fasta)
        return new_fasta
    else:
        clean_cap3(fasta)
        return fasta

def clean_cap3(fasta):
    #Clean extra files
    for ext in ['.cap.contigs.links', '.cap.contigs.qual', '.cap.info', '.cap.ace', '.cap.singlets', '.cap.contigs']:
        if os.path.isfile('%s%s' % (fasta, ext)):
            os.remove('%s%s' % (fasta, ext))


def join_cap3_contigs(fasta):
    sing_names, sing_seqs = read_fasta_lists_simple_names('%s.cap.singlets' % fasta)
    joined_names, joined_seqs = read_fasta_lists_simple_names('%s.cap.contigs' % fasta)
    all_names=sing_names+joined_names
    all_seqs=sing_seqs+joined_seqs
    cap3_cons_name='%s/cap3_%s' % ('/'.join(fasta.split('/')[:-1]), fasta.split('/')[-1])
    write_fasta(all_names, all_seqs, cap3_cons_name)
    return cap3_cons_name

def cap3_improved(fasta):
    old_names, old_seqs = read_fasta_lists_simple_names(fasta)
    new_sing_names, new_sing_seqs = read_fasta_lists_simple_names('%s.cap.singlets' % fasta)
    if len(new_sing_names) >= len(old_names): return False
    else: return True

# will cut fasta name off at the first whitespace
def read_fasta_lists_simple_names(file):
    fin = open(file, 'r')
    count=0

    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:].split()[0])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)

    return names, seqs

def run_blast(opts):
    #Create symbolic link in the cwd to the contigs file
    query='%s/ray_contigs.fasta' % (os.getcwd())
    if not os.path.exists(query):
        os.symlink('%s/Contigs.fasta' % (opts.rayOut), query)
    else: print 'Problem! There is already a file with the name "ray_contigs.fasta"'
    
    if opts.ns and opts.ps: cmd='%s -q %s -b %s -n %d --ns %s --ps %s  --temp %s --task %s --evalue %s --goodHit %s --numHits %d --numHsps %d' % (opts.blast, query, opts.blastType, opts.procs, opts.ns, opts.ps, opts.temp, opts.task, opts.evalue, opts.evalue, opts.numHits, opts.numHsps)
    elif opts.ns: cmd='%s -q %s -b %s -n %d --ns %s  --temp %s --task %s --evalue %s --goodHit %s --numHits %d --numHsps %d' % (opts.blast, query, opts.blastType, opts.procs, opts.ns, opts.temp, opts.task, opts.evalue, opts.evalue, opts.numHits, opts.numHsps)
    else: cmd='%s -q %s -b %s -n %d --ps %s  --temp %s --evalue %s --goodHit %s --numHits %d --numHsps %d' % (opts.blast, query, opts.blastType, opts.procs, opts.ps, opts.temp, opts.evalue, opts.evalue, opts.numHits, opts.numHsps)

    #Add the dontIndex flag, if thrown
    if opts.dontIndex: cmd+=' --dontIndex'

    #Add in redirection of stdin and stdout
    cmd+=' &>>blast_stderr'
    
    print cmd
    blastit=Popen(cmd, shell=True)
    blastit.wait()
    
    names_interest_dict = parse_blast(query, opts)
    if len(names_interest_dict)>0:
        fasta_interest = 'contigs_interest.fasta'
        subset_fasta(names_interest_dict, query, fasta_interest, opts)
        return fasta_interest
    else: return None

def parse_blast(query, opts):
    names_interest={}
    types_tasks_tups = get_types_run(opts)
    for info in types_tasks_tups: 
        if info[0]=='blastn' or info[0]=='tblastx': subject = opts.ns
        else: subject = opts.ps
        fin=open('%s_%s_%s_%s_parsed.txt' % (query, info[0], info[1][:2], subject.split('/')[-1]), 'r')
        for line in fin:
            cols=line.strip().split('\t')
            names_interest[cols[0]]=''
            query='%s_%s_no_good_hits.fasta' % (info[0], info[1])
    return names_interest

def get_types_run(opts):
    base_types=opts.blastType.split(',')
    types_tasks=[]
    for ty in base_types: 
        if ty=='blastn':
            base_tasks=opts.task.split(',')
            for ta in base_tasks: types_tasks.append((ty, ta))
        else: types_tasks.append((ty, ''))
    return types_tasks

def run_ray2(base_input, opts):
    cmd='mpiexec -np %d %s -disable-scaffolder -k %d -o %s -p %s %s -s %s -s %s &>>ray_stderr' % (opts.procs, opts.ray, opts.kmer, opts.rayOut, '%s_1.fastq' % base_input, '%s_2.fastq' % base_input, '%s_1_singletons.fastq' % base_input, '%s_2_singletons.fastq' % base_input)
    print cmd
    cutit=Popen(cmd, shell=True)
    cutit.wait()

def rev_comp(seq):
    dna = Seq(seq, generic_dna)
    return str(dna.reverse_complement())

#Use this to delete some files, given that those files exist    
def delete_files(*done_files):
    for to_remove in done_files:
        if os.path.isfile(to_remove):
            os.remove(to_remove)

def run_cutadapt(fastq, ill_adapt, opts):
    if fastq.endswith('.gz'): out_name = '%s_cutadapt.fastq' % '.'.join(fastq.split('.')[:-2])
    else: out_name = '%s_cutadapt.fastq' % '.'.join(fastq.split('.')[:-1])
    cmd='%s -a %s -g %s -a %s -o %s -m %d --match-read-wildcards %s' % (opts.ca, ill_adapt, opts.sispa, opts.sispa_revcomp, out_name, opts.minLength, fastq)
    print cmd
    cutit=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    cutit.wait()
    for line in cutit.stdout:
        print line.rstrip()
    return out_name
    
def run_prinseq(r1, r2, opts):
    out_name = '%s_good' % '.'.join(r1.split('.')[:-1])
    if not opts.NOderep: cmd='%s -derep %s -out_good %s -out_bad null -min_len %d -lc_method %s -lc_threshold %s -trim_right %d -min_qual_mean %d -fastq %s -fastq2 %s' % (opts.pseq, opts.derep, out_name, opts.minLength, opts.lcMethod, opts.lcThresh, opts.trimRight, opts.meanQual, r1, r2)
    else: cmd='%s -out_good %s -out_bad null -min_len %d -lc_method %s -lc_threshold %s -trim_right %d -min_qual_mean %d -fastq %s -fastq2 %s' % (opts.pseq, out_name, opts.minLength, opts.lcMethod, opts.lcThresh, opts.trimRight, opts.meanQual, r1, r2)
    print cmd
    ps=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    ps.wait()
    for line in ps.stderr:
        print line.rstrip()
        
    return out_name
    
#For sifting through fastqs and just keeping what is still paired

def paired_sort(read1, read2):
    names1 = get_names(read1)
    names2 = get_names(read2)
    paired = set(names1) & set(names2)
	
    del names1
    del names2
	
    pair1_file = write_new_file(read1, paired)
    pair2_file = write_new_file(read2, paired)
    
    return pair1_file, pair2_file

def get_names(file):
    fin = open(file, 'r')	
    names=[]
    linenum=0
	
    for line in fin:
        linenum+=1		
        #First name line
        if linenum%4==1:
            names.append(line.strip().split()[1])
    fin.close()
    return names

def write_new_file(fastq, paired_names):

    fin = open(fastq, 'r')
    fout_pair = open(fastq.split('.')[0] + '_paired.fastq', 'w')
    linenum=0
    is_paired=0
	
    for line in fin:
        linenum+=1
        #First name line
        if linenum%4==1:
            name=line.strip().split()[1]
            if name in paired_names:
                is_paired=1
                fout_pair.write(line)	
            else:
                is_paired=0
        #Other lines
        else:
            if is_paired: fout_pair.write(line)
    fin.close()
    fout_pair.close()
    return fastq.split('.')[0] + '_paired.fastq'

###------------->>>

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
	
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
	
    return names, seqs

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()

def subset_fasta(subset, fasta, new_name, opts):
    names, seqs=read_fasta_lists(fasta)
    nametrim=[name.rstrip() for name in names]
    subnames=[]
    subseqs=[]
	
    for index in range(len(names)):
        if nametrim[index] in subset and len(seqs[index])>=opts.minContigLen:
            subnames.append(names[index])
            subseqs.append(seqs[index])

    write_fasta(subnames, subseqs, new_name)

            
###------------>>>

if __name__ == "__main__":
    main()	
