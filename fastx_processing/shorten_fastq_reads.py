#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fastq, can\'t be gzipped [None,REQD]')
  p.add_option('-l', '--length', type='int', default= 50, help='Length to keep [50]')
  p.add_option('-o', '--output', help='Output fastq [None, REQD]')
  
  opts, args = p.parse_args()

  with open(opts.input, "r") as fin:
    with open(opts.output, "w") as fout:
      for record in FastqGeneralIterator(fin):
        fout.write("@%s\n%s\n+\n%s\n" % (record[0], record[1][:opts.length], record[2][:opts.length]))

if __name__ == '__main__':
  main()
