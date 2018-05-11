#!/usr/local/bin/python2.7

#This script reads a file that contains info for a number of subsequences that you want to output to a new file
#The new seqs are output to stdout
#Each line of the file specified after the command line should include: fasta_file, seq_name, base_start, base_stop

import sys
from Bio.Seq import Seq

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

def read_fasta_dict(file):
    names, seqs = read_fasta_lists(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict

def get_seq(file, seq_name, start, stop, strand):
    seq_dict = read_fasta_dict(file)
    seq_of_interest = seq_dict[seq_name][int(start)-1:int(stop)]
    if strand == "-":
      forward_seq = Seq(seq_of_interest)
      seq_of_interest = forward_seq.reverse_complement()
#    print '>%s_%d_%d %s \n%s' % (seq_name, int(start), int(stop), file, seq_of_interest)
    return seq_of_interest

###---------->>>

if __name__=='__main__':
    if len(sys.argv)==2:
        fin = open(sys.argv[1], 'r')
        for line in fin:
            cols = line.strip().split('\t')
            file_name = cols[0]
            ref_name = cols[1]
            seq_name = cols[2]
            seq = ""
            for i in range(3,len(cols)-2,3):
              seq += get_seq(file_name, ref_name, cols[i], cols[i+1], cols[i+2])
            print '>%s %s\n%s' % (seq_name, file_name, seq)

    else: print 'Each line of the file specified after the command line should include: fasta_file, ref_name, seq_name, base_start, base_stop, orientation. You can repeat the last 3 columns as necessary'
