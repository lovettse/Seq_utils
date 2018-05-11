#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
from heapq import nlargest
import random
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fastq [None,REQD]')
  p.add_option('-n', '--number', type='int', default= 50000, help='Numnber of reads to keep [50000]')
  p.add_option('-o', '--output', help='Output fastq [None,REQD]')
  
  opts, args = p.parse_args()

  with open(opts.input, "r") as fin:
    with open(opts.output, "w") as fout:
      # TBH, I don't understand how this works. Check here for an explanation:
      # stackoverflow.com/questions/12581437/python-random-sample-with-a-generator-iterable-iterator
      for record in (x for _, x in nlargest(opts.number, ((random.random(),x) for x in FastqGeneralIterator(fin)))):
        fout.write("@%s\n%s\n+\n%s\n" % (record[0], record[1], record[2]))

if __name__ == '__main__':
  main()