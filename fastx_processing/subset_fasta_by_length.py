#!/usr/bin/env python

from Bio.SeqIO.FastaIO import SimpleFastaParser
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fasta [None,REQD]')
  p.add_option('--min', type="int", default=1000, help="Minimum size of record to keep [1000]")
  p.add_option('--max', type="int", help="Maximum size of record to keep [None]")
  p.add_option('-o', '--output', help='Output fasta [None,REQD]')
  
  opts, args = p.parse_args()

  with open(opts.input, "r") as fin:
    with open(opts.output, "w") as fout:
      for record in SimpleFastaParser(fin):
        if opts.max:
          if len(record[1]) >= opts.min and len(record[1]) <= opts.max:
            fout.write(">%s\n%s\n" % (record[0], record[1]))
        else:
          if len(record[1]) >= opts.min:
            fout.write(">%s\n%s\n" % (record[0], record[1]))


if __name__ == '__main__':
  main()
