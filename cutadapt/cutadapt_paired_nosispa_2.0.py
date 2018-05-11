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
#    p.add_option('--ps', default='prinseq-lite.pl', help='How to call prinseq. [prinseq-lite.pl]')
    
    #Settings
    p.add_option('--i1', default='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC', help="Illumina adaptor to be trimmed from the 3' end of R1. [GATCGGAAGAGCACACGTCTGAACTCCAGTCAC]")
    p.add_option('--i2', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT', help='Reverse complement of the Illumina adaptor that will be trimmed from the 3; end of R2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')
    p.add_option('-l', '--minLength', default=70, type='int', help='Minimum length to keep (after primer clipping) [70]')
    p.add_option('-q', '--quiet', action="store_true",dest="quiet",default=False,help='Quiet mode [False]')
#    p.add_option('-q', '--meanQual', default=20, type='int', help='Minimum mean quality for a read to be kept [20]')
    
    opts, args = p.parse_args()
    
    #Run cutadapt
    one_cut, two_cut = run_both_cutadapts(opts)


    #Make new fastqs with just reads that are still paired. Delete intermediate files
    one_cut_paired, two_cut_paired = paired_sort(one_cut, two_cut)
    delete_files(one_cut, two_cut)	

#End of main()

def run_both_cutadapts(opts):
    if opts.one.endswith('.gz'): out1_name = '%s_cutadapt.fastq' % '.'.join(opts.one.split('.')[:-2])
    else: out1_name = '%s_cutadapt.fastq' % '.'.join(opts.one.split('.')[:-1])
    cmd='%s -a %s -o %s -m %d --match-read-wildcards %s' % (opts.ca, opts.i1, out1_name, opts.minLength, opts.one)
    if not opts.quiet: print cmd
    cutit1=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)

    if opts.two.endswith('.gz'): out2_name = '%s_cutadapt.fastq' % '.'.join(opts.two.split('.')[:-2])
    else: out2_name = '%s_cutadapt.fastq' % '.'.join(opts.two.split('.')[:-1])
    cmd='%s -a %s -o %s -m %d --match-read-wildcards %s' % (opts.ca, opts.i2, out2_name, opts.minLength, opts.two)
    if not opts.quiet: print cmd
    cutit2=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)    

    cutit1.wait()
    cutit2.wait()

    if not opts.quiet:
        for line in cutit1.stdout:
            print line.rstrip()
        for line in cutit2.stdout:
            print line.rstrip()
    return out1_name, out2_name


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
    
def batch_filter(opts, cut_files):
    outfiles=[]
    for file in cut_files:
        out_name = '%s_good' % '.'.join(file.split('.')[:-1])
        outfiles.append(out_name)
        cmd='prinseq-lite.pl -derep 14 -out_good %s -out_bad null -min_len %d -lc_method dust -lc_threshold 3 -trim_right 20 -min_qual_mean %d -fastq %s' % (out_name, opts.minLength, opts.meanQual, file)
        ps=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        ps.wait()
        for line in ps.stderr:
            print line.rstrip()
        
    return outfiles
    
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
    return fout_pair

###------------->>>


            
###------------>>>

if __name__ == "__main__":
    main()	
