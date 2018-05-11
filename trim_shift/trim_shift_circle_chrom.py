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
                                           Each line should have either 3 or 5 tab-dleimited columns.
                                                - The first 3 should always be: 
                                                        contig_name
                                                        position_of_new_start_base (pre-trim)
                                                        a number idicating whether the sequence should be reverse complemented (0=no, 1=yes) 
                                                - If present, the fourth and fifth should be:
                                                        the first base to keep from the beginning of the sequence
                                                        the last base to keep from the end of the sequence
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
    
    for name, seq in input_dict.iteritems():
        #If some adjustments need to be made
        if name in to_change_dict:
            seq_copy=seq[:]
            front_ind=to_change_dict[name][2]
            end_ind=to_change_dict[name][3]
            new_start=to_change_dict[name][0]
            #trims ends, if specified. Also adjusts the new start based on the number of bases trmmed from the begining.
            if len(to_change_dict[name])==4:
                seq_copy=seq_copy[front_ind-1:end_ind]
                new_start=new_start-front_ind+1
                #print seq_copy, new_start
            #Reverse complement, if needed
            if to_change_dict[name][1]:
                seq_copy=rev_comp(seq_copy)
                new_start=len(seq_copy)-new_start+1
            #Shift chromosomal break point
            seq_copy=rebreak(seq_copy, new_start)

            input_dict[name]=seq_copy
    new_names=input_dict.keys()
    new_seqs=[input_dict[x] for x in new_names]
    write_fasta(new_names, new_seqs, opts.out)
                
def rebreak(seq, new_first):
    new_seq=seq[new_first-1:] + seq[:new_first-1]
    return new_seq

# This should look fine regardless of whether there are 3 columns (no trimming) or five columns (trimming)
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
