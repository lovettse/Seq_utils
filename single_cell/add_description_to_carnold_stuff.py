#!/usr/bin/env python

from Bio import SeqIO
from collections import Counter
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input file (whatever format this is) [None, REQD]')
  p.add_option('-f', '--fasta', help='Input fasta [None,REQD]')
  p.add_option('-o', '--output', help='Output file [None,REQD]')
  
  opts, args = p.parse_args()

  # key = BAT... id, value = LOC... id
  fasta_dict = {}
  out_lines = []

  with open(opts.fasta, "r") as fin:
    for line in fin:
      if line.startswith(">"):
        cols = line[1:].strip().split("\t")
        fasta_dict[cols[0]] = cols[1]

  with open(opts.input, "r") as fin:
    header = fin.readline().replace('"', '')
    header = header.split(" ")[0] + "\tproduct\t" + "\t".join(header.split(" ")[1:])

    for line in fin:
      cols = line.strip().replace('"', '').split(" ")
      new_line = "\t".join(cols[:2]) + "\t%s\t" % fasta_dict[cols[1]] + "\t".join(cols[2:]) + "\n"
      out_lines.append(new_line)

  with open(opts.output, "w") as fout:
    fout.write("%s" % header)
    for line in out_lines:
      fout.write(line)
    

if __name__ == '__main__':
  main()