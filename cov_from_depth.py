#!/usr/bin/python

from __future__ import division
import optparse

def main():


  usage = '%prog [options] > report.txt'
  p = optparse.OptionParser(usage=usage)

  p.add_option('-i', '--input', help='Input depth report, must be run with "samtools depth -a" [None, REQD]')

  opts, args = p.parse_args()
  
  total_bases = 0
  ref_len = 0
  zero_cov = 0
  
  with open(opts.input, "r") as fin:
    for line in fin:
      this_depth = int(line.strip().split("\t")[2])
      total_bases += this_depth
      if this_depth == 0:
        zero_cov += 1
     
      ref_len += 1

#      if ref_len % 1000000 == 0:
#        print "%s lines processed" % (ref_len)

  print "\n"
  print "Length of reference: %s" % ref_len
  print "Total bases: %s" % total_bases
  print "Average depth: %s" % (total_bases/ref_len)
  print "Zero coverage bases: %s" % zero_cov
  print "Percent zero coverage: %s" % (zero_cov/ref_len)

if __name__=="__main__":
  main()
