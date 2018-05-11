#!/usr/bin/env python

# This script checks for identical sequences in an alignment
# It does not calculate values for pairwise identity

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
  iden_dict = {}
  final_dict = {}
  
  with open(opts.input, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      fasta_dict[record.id] = str(record.seq).replace("n", "X").replace("N","X").replace("-","X")
    align_len = len(record.seq)

  for header1 in fasta_dict:
    iden_dict[header1] = []
    for header2 in fasta_dict:
      if not header1 == header2 and fasta_dict[header1] == fasta_dict[header2]:
        if not header1 in iden_dict:
          iden_dict[header1] = []
        iden_dict[header1].append(header2)
    iden_dict[header1] = sorted(iden_dict[header1])

  for header1 in iden_dict:
    if iden_dict[header1]:
      for header2 in final_dict:
        test1 = iden_dict[header1][:]
        test1.append(header1)
        test1 = sorted(test1)
        test2 = final_dict[header2][:]
        test2.append(header2)
        test2 = sorted(test2)
        if test1 == test2:
          break
      else:
        final_dict[header1] = iden_dict[header1]

  with open(opts.output, "w") as fout:
    for header in final_dict:
      fout.write("%s\t%s\n" % (header, ",".join(final_dict[header])))

if __name__ == '__main__':
  main()
