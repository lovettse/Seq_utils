#!/usr/bin/env python

import optparse

def main():
  usage = '%prog [options]'
  p = optparse.OptionParser()
  p.add_option('-i', '--input', help="Input table [None, REQD]")
  p.add_option('-o', '--output', help="Output table [None, REQD]")

  opts, args = p.parse_args()

  with open(opts.input,"r") as fin:
    lis = [x.strip().split("\t") for x in fin]

  with open(opts.output, "w") as fout:
    for x in zip(*lis):
      for y in x:
        fout.write("%s\t" % y)
      fout.write("\n")

if __name__=="__main__":
  main()
