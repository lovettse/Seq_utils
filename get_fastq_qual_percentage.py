#!/usr/bin/env python

from __future__ import division
from Bio import SeqIO
import optparse, gzip

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input fastq. Can be gzipped if the ".gz" extension is present. [None, REQD]')
  p.add_option('-q', '--qual', type="int", default=30, help='Quality level [30]')
  p.add_option('-o', '--output', default="/dev/stdout", help="Output file. Will print to screen if not specified.\nFile will be appended to, so erase is it you want to start fresh [/dev/stdout]")

  opts, args = p.parse_args()

  total_records = 0
  pass_filter = 0

  if opts.input.endswith(".gz"):
    with gzip.open(opts.input, "rt") as fin:
      for record in SeqIO.parse(fin, "fastq"):
        total_records += 1
        scores = record.letter_annotations["phred_quality"]
        if sum(scores)/len(scores) > opts.qual:
          pass_filter += 1
  else:
    for record in SeqIO.parse(opts.input, "fastq"):
      total_records += 1
      scores = record.letter_annotations["phred_quality"]
      if sum(scores)/len(scores) > opts.qual:
        pass_filter += 1

  with open(opts.output, "a") as fout:
    fout.write("%s\t%0.4f\n" % (opts.input.split("/")[-1],pass_filter/total_records*100))

if __name__=="__main__":
  main()