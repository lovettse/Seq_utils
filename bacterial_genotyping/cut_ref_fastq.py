#!/usr/bin/env python

####################################################################################################
# By Jason Ladner
# Takes a 'reference' fasta and splits it into 'reads' of a specified length, tiling across the genome
# Default is to have a step size of 1 so that there is one read per starting base, until the end when the reads would be too short
# The names of the new reads are based on the reference sequence name and the starting base position 
# The reads are output in fastq format with base qualities all set to 32
###################################################################################################

from __future__ import division
import sys, optparse, os, math

def main():

    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)

    p.add_option('-f', '--fasta', help='Input "reference" fasta file to be cut up into reads [None, REQ]')
    p.add_option('-w', '--win', type='int', default=70, help='Size for the reads [70]')
    p.add_option('-s', '--slide', type='int', default=1, help='Offset in bases between reads [1]')
    p.add_option('-o', '--out', help="Name for output fastq file [None, REQ]")
    
    #All extra arguments should be databases formatted for blast
    opts, args = p.parse_args()    

    cut_fasta(opts)
#------------------------->>>

def cut_fasta(opts):
    #Read in fasta file
    names, seqs = read_fasta_lists(opts.fasta)
    
    #Open output file for writing
    fout=open(opts.out, 'w')
    
    #Step through each ref seq and build the subset fastq reads
    for index, each in enumerate(names):
        j=0
        while j<len(seqs[index])-opts.win+1:
            fout.write('@%s_%d\n%s\n+\n%s\n' % (each, j+1, seqs[index][j:j+opts.win], 'A'*opts.win))
            j+=opts.slide
    
    fout.close()
    
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


def calc_div_size(lines, divs):
    num_seqs = lines/4
    sub_size = math.ceil(num_seqs/divs)
    return int(sub_size)

def count_lines(file):
    fin = open(file, 'r')
    line_count=0
    for line in fin: line_count+=1
    fin.close()
    return line_count 

def make_dir_if_necc(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def divide_fastq(o, div_size):
    fin = open(o.inputFastq, 'r')
    basename_out = '.'.join(o.inputFastq.split('.')[:-1])
    total_seq_count=0
    while 1:
        fout = open('%s/%s_%d_%d.fastq' % (o.outDir, basename_out, total_seq_count+1, total_seq_count+div_size), 'w')
        seq_count=0
        while seq_count<div_size:
            total_seq_count+=1
            seq_count+=1
            
            seq = []
            #Read Next Sequence
            for i in range(4):
                seq.append(fin.readline())
            
            #Check to make sure you are not at the end of the file
            if seq[0]!='': 
                fout.writelines(seq)
            
            else:
                #Close all open files
                fout.close()
                fin.close()
                return
            
        fout.close()
        
###------------->>>

if __name__ == "__main__":
    main()
