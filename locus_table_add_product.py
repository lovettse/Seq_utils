#!/usr/bin/env python

'''
By Sean Lovett
9/27/2016
'''

from __future__ import print_function
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input locus table [None, REQD]')
  p.add_option('-p', '--prods', help='Input list of locus tags and product names, tab delimited [None, REQD]')
  p.add_option('-o', '--output', help='Output ortholog table [None,REQD]')
  
  opts, args = p.parse_args()
  
  prod_dict = {}
  
  with open(opts.prods, "r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      prod_dict[cols[0]] = cols[1]

  with open(opts.input, "r") as fin:
    with open(opts.output, "w") as fout:
      for line in fin:
        outs = []
        cols = line.strip().split("\t")
        for col in cols:
          genes = col.split(',')
          for gene in genes:
            if gene in prod_dict:
              outs.append("%s %s" % (gene,prod_dict[gene]))
        for out in outs:
          fout.write("%s;\t" % out)
        fout.write("\n")


if __name__ == "__main__":
  main()