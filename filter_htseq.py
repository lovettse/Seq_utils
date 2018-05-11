#!/usr/bin/env python

'''
By Sean Lovett
11/3/2016
'''

from __future__ import print_function
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input htseq file [None, REQD]')
  p.add_option('-l', '--list', help='List of locus tags, can have tab-delimited columns [None, REQD]')
  p.add_option('-o', '--output', help='Output htseq file, will only contain locus tags that were in the list [None,REQD]')

  opts, args = p.parse_args()

  in_dict = {}
  tag_list = []

  with open(opts.input,"r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      in_dict[cols[0]] = cols[1]

  with open(opts.list, "r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      for col in cols:
        tag_list.append(col)

  with open(opts.output, "w") as fout:
    for tag in sorted(tag_list):
      if tag in in_dict:
        fout.write("%s\t%s\n" % (tag, in_dict[tag]))

if __name__=="__main__":
  main()