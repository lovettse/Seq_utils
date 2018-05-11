#!/usr/bin/env python

import optparse

def main():
  usage = '%prog [options]'
  p = optparse.OptionParser()
  p.add_option('-i', '--input', help="Input ortho table [None, REQD]")
  p.add_option('-o', '--output', help="Output file [None, REQD]")

  opts, args = p.parse_args()

  lines = []

  with open(opts.input, "r") as fin:
    header = fin.readline().strip().split("\t")[1:]
    for line in fin:
      cols = line.strip().split("\t")[1:]
      count = 0
      for col in cols:
        if col and not "," in col:
          count += 1
      if count == len(header):
        lines.append("\t".join(cols))

  with open(opts.output, "w") as fout:
    fout.write("%s\n" % "\t".join(header))
    for line in lines:
      fout.write("%s\n" % line)

if __name__=="__main__":
  main()
