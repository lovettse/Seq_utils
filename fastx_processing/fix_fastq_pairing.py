#!/usr/bin/python

from __future__ import division
import optparse
from Bio import SeqIO

def main():


  usage = '%prog [options]'
  p = optparse.OptionParser()

  p.add_option('-1', '--read1', help='Input R1 [None, REQD]')
  p.add_option('-2', '--read2', help='Input R2 [None, REQD]')

  opts, args = p.parse_args()
  
  paired_sort(opts.read1, opts.read2)    

    
  


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


def paired_sort(read1, read2):
    names1 = get_names(read1)
    names2 = get_names(read2)
    paired = set(names1) & set(names2)

    del names1
    del names2

    pair1_file = write_new_file(read1, paired)
    pair2_file = write_new_file(read2, paired)

    return pair1_file, pair2_file

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


if __name__=="__main__":
  main()
