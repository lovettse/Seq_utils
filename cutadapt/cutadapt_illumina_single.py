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
    
    #How to call programs
    p.add_option('--ca', default='cutadapt1.2', help='How to call cutadapt. Need to use version >=1.2 [cutadapt1.2]')
    p.add_option('--pseq', default='prinseq-lite_0.20.3.pl', help='How to call prinseq. [prinseq-lite_0.2.2.pl]')
    p.add_option('--ray', default='Ray2.2', help='How to call Ray. [Ray2.2]')
    p.add_option('--cap', default='cap3', help='How to call cap3. [cap3]')
    p.add_option('--bb', default='bowtie_buildout_latest_wQualchecks.py', help='How to call bowtie_buildout. [bowtie_buildout_latest_wQualchecks.py]')
    p.add_option('--blast', default='split_blasts_multi_faster.py', help='How to call blast function. [split_blasts_multi_faster.py]')
    
    #Cutadapt Settings
    p.add_option('--i1', default='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', help="Illumina adaptor to be trimmed from the 3' end of R1. [GATCGGAAGAGCACACGTCTGAACTCCAGTCAC]")
    p.add_option('--i2', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT', help='Reverse complement of the Illumina adaptor that will be trimmed from the 3; end of R2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')
    
    #Prinseq Settings
    p.add_option('-l', '--minLength', default=50, type='int', help='Minimum length to keep (after primer clipping) [50]')
    p.add_option('-q', '--meanQual', default=20, type='int', help='Minimum mean quality for a read to be kept [20]')
    p.add_option('--derep', default='14', help='derep setting for prinseq. Use the --NOderep flag to turn replicate removal off [14]')
    p.add_option('--NOderep', default=False, action='store_true',  help='Use this flag if you do not want to remove duplicates [False]')
    p.add_option('--lcMethod', default='dust',  help='Tells prinseq which method to use to get rid of low complexity sequences [dust]')
    p.add_option('--lcThresh', default='3',  help='Low complexity threshold [3]')
    p.add_option('--trimRight', default=20, type='int',  help='All bases below this quality will be trimmed from the right of the sequences. Must be an integer [20]')


    
    opts, args = p.parse_args()
    
    #Run cutadapt
    one_cut = run_cutadapt(opts.one, opts.i1, opts)
    #two_cut = run_cutadapt(opts.two, opts.i2, opts)

    #Make new fastqs with just reads that are still paired. Delete intermediate files
    #one_cut_paired, two_cut_paired = paired_sort(one_cut, two_cut)
    #delete_files(one_cut, two_cut)	
        
    base_good_name = run_prinseq_single(one_cut, opts)

#-------------------------------End of main()------------------------


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
    cmd='%s -a %s -o %s -m %d --match-read-wildcards %s' % (opts.ca, ill_adapt, out_name, opts.minLength, fastq)
    print cmd
    cutit=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    cutit.wait()
    for line in cutit.stdout:
        print line.rstrip()
    return out_name
    
def run_prinseq_single(r1, opts):
    out_name = '%s_good' % '.'.join(r1.split('.')[:-1])
    if not opts.NOderep: cmd='%s -derep %s -out_good %s -out_bad null -min_len %d -lc_method %s -lc_threshold %s -trim_right %d -min_qual_mean %d -fastq %s' % (opts.pseq, opts.derep, out_name, opts.minLength, opts.lcMethod, opts.lcThresh, opts.trimRight, opts.meanQual, r1)
    else: cmd='%s -out_good %s -out_bad null -min_len %d -lc_method %s -lc_threshold %s -trim_right %d -min_qual_mean %d -fastq %s' % (opts.pseq, out_name, opts.minLength, opts.lcMethod, opts.lcThresh, opts.trimRight, opts.meanQual, r1)
    print cmd
    ps=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    ps.wait()
    for line in ps.stderr:
        print line.rstrip()
        
    return out_name
    
#For sifting through fastqs and just keeping what is still paired


def get_names(file):
    fin = open(file, 'r')	
    names=[]
    linenum=0
	
    for line in fin:
        linenum+=1		
        #First name line
        if linenum%4==1:
            names.append(line.strip().split()[0])
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
            name=line.strip().split()[0]
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
