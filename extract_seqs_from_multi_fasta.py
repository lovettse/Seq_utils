#!/usr/bin/env python

from __future__ import division
from Bio import SeqIO
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-f', '--fasta', help='Fasta to get seqs from [None,REQD]')
  p.add_option('-l', '--list', help='List of seqs to retrieve [None, REQD]')
  p.add_option('-o', '--output', help='Output fasta [None, REQD]')

  opts, args = p.parse_args()

  seq_names = []
  seq_dict = {}

  with open(opts.list, "r") as fin:
    for line in fin:
      seq_names.append(line.strip())

  with open(opts.fasta, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      if record.id in seq_names:
        seq_dict[record.id] = str(record.seq)

  with open(opts.output, "w") as fout:
    for seq in seq_dict:
      fout.write(">%s\n%s\n" % (seq, seq_dict[seq]))

if __name__=="__main__":
  main()
