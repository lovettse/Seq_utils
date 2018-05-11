#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fastq [None,REQD]')
  p.add_option('-s', '--size', type="int", default=1000, help="Minimum size of read to keep [1000]")
  p.add_option('-o', '--output', help='Output fastq [None,REQD]')
  
  opts, args = p.parse_args()

  with open(opts.input, "r") as fin:
    with open(opts.output, "w") as fout:
      for record in FastqGeneralIterator(fin):
        if len(record[1]) >= opts.size:
          fout.write("@%s\n%s\n+\n%s\n" % (record[0], record[1], record[2]))


if __name__ == '__main__':
  main()