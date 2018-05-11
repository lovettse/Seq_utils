#!/usr/bin/env python

from Bio import SeqIO
from collections import Counter
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fasta [None,REQD]')
  
  opts, args = p.parse_args()

  with open(opts.input, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      seq = Counter(record.seq.upper())
      print "%s\t%s" % (record.id, seq['N'] + seq['-'])

if __name__ == '__main__':
  main()