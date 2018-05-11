#!/usr/bin/env python

'''
By Sean Lovett
9/6/2016
'''

from __future__ import print_function
from subprocess import Popen, PIPE
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input ffn. Duplicate headers allowed. [None, REQD]')
  p.add_option('-o', '--output', help='Output ffn. Will only contain the longest sequence with each header. Shorter sequences with the same header removed. [None,REQD]')
  
  opts, args = p.parse_args()
  
  fasta_dict = read_fasta_dict(opts.input)
  write_fasta_from_dict(fasta_dict, opts.output)
  
def read_fasta_dict(file):
  names, seqs = read_fasta_lists(file)
  fasta_dict = dict(zip(names, seqs))
  return fasta_dict

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
      names.append(line[1:].strip().split("|")[1])
      if count>1:
        seqs.append(seq.strip())
      seq=''
    else:
      seq += line
  seqs.append(seq)

  return names, seqs

def write_fasta_from_dict(seq_dict, out_file):
  with open(out_file,"w") as fout:
    for key, seq in seq_dict.items():
      fout.write(">%s\n%s\n" % (key, seq))
  
  
if __name__ == "__main__":
  main()