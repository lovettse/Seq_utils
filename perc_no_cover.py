#!/usr/bin/env python

from __future__ import division
from Bio import SeqIO
import optparse
import re

def main():
  #To parse command line
  usage = "usage: %prog -i input.txt -r ref.fasta -g ref.gff [options] > commandlog.txt"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fasta [None,REQD]')
  p.add_option('-o', '--output', help='Output fasta [None,REQD]')
  
  opts, args = p.parse_args()
  
  fasta_dict = {}
  out_dict = {}
  
  with open(opts.input, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      fasta_dict[record.id] = str(record.seq).replace("n", "X").replace("N","X").replace("-","X")
    align_len = len(record.seq)

  for header in fasta_dict:
    blanks = 0
    for i in range(len(fasta_dict[header])):
      if fasta_dict[header][i] == "X":
        blanks += 1
    out_dict[header] = blanks
    
  with open(opts.output, "w") as fout:
    for header in sorted(out_dict):
      fout.write("%s\t%0.8f\n" % (header, out_dict[header]/align_len)) 

    


if __name__ == '__main__':
  main()