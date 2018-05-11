#!/usr/bin/env python

from __future__ import division
import optparse, os
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--input', help='Input fasta with contig sequence(s) [None, REQD]')
    p.add_option('-o', '--out', default="out_circ_contigs.fasta", help='Name for the output file [out_circ_contigs.fasta]')
    p.add_option('-b', '--breaks', help='''File with new breakpoint info for each contig that needs to be reoriented.
                                           This file should contain one line for each sequence to be altered.
                                           Each line should be 5 columns.
                                                contig_name
                                                start of backbone sequence
                                                end of backbone sequence
                                                a number idicating whether the sequence should be reverse complemented (0=no, 1=yes) 
                                                amount of circular redundancy
                                                [None, REQ]''')

    opts, args = p.parse_args()
    
    rearrange_chroms(opts)

###--------------end of main()------------------

def rev_comp(seq):
    dna = Seq(seq, generic_dna)
    return str(dna.reverse_complement())

def rearrange_chroms(opts):
    input_dict = read_fasta_dict(opts.input)
    to_change_dict=get_break_info(opts.breaks)
    
    new_names = []
    
    for name, seq in input_dict.iteritems():
        #If some adjustments need to be made
        if name in to_change_dict:
            new_names.append(name)
            seq_copy=seq[:]
            start_backbone = to_change_dict[name][0]
            end_backbone = to_change_dict[name][1]
            redundant_count = to_change_dict[name][3]
            to_trim = int(redundant_count/2)
            
            #trims ends, if specified. Also adjusts the new start based on the number of bases trmmed from the begining.
            print ("I should get %s to %s and %s to %s" % (to_trim-1, start_backbone, end_backbone, len(seq_copy)-to_trim))
            seq_first_piece = seq_copy[to_trim:start_backbone]
            seq_second_piece = seq_copy[end_backbone:len(seq_copy)-to_trim]

            seq_copy = seq_second_piece + seq_first_piece
            
            #Reverse complement, if needed
            if to_change_dict[name][2]:
                seq_copy=rev_comp(seq_copy)

            input_dict[name]=seq_copy
    #new_names=input_dict.keys()
    new_seqs=[input_dict[x] for x in new_names]
    write_fasta(new_names, new_seqs, opts.out)
                
def rebreak(seq, new_first):
    new_seq=seq[new_first-1:] + seq[:new_first-1]
    return new_seq


def get_break_info(bfile):
    to_change_dict={}
    fin=open(bfile, 'r')
    for line in fin:
        cols=line.strip().split('\t')
        to_change_dict[cols[0]]=[int(x) for x in cols[1:]]
    return to_change_dict


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

def read_fasta_dict(file):
    names, seqs = read_fasta_lists(file)
    fasta_dict = dict(zip(names, seqs))
    return fasta_dict


###---------->>>

if __name__ == '__main__':
    main()
