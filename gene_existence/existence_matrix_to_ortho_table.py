#!/usr/bin/env python

'''
By Sean Lovett
9/12/2016
'''

from __future__ import print_function
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input existence matrix. Output may need trimming if not all genes are present in all species [None, REQD]')
  p.add_option('-o', '--output', default="ortho_table.txt", help='Output ortholog table [ortho_table.txt]')
  p.add_option('--exclude', default=False, action="store_true", help="Exclude ortholog groups without a member in every genome [False]")
  
  opts, args = p.parse_args()
  
  with open(opts.input,"r") as fin:
    with open(opts.output,"w") as fout:
      header = fin.readline().strip().split("\t")[1:-1]
      fout.write("\t%s\n" % "\t".join(header))
      for in_line in fin:
        not_all_present = False
        out_line = []
        in_cols = in_line.strip().split("\t")[1:-1]
        name = in_line.strip().split("\t")[0]
        out_line.append(name)
        for i in range(len(in_cols)):
          if in_cols[i] == "1":
            out_line.append("_".join((header[i],name)))
          else:
            not_all_present = True
            out_line.append("NA")
        if opts.exclude and not_all_present:
          pass
        else:
          fout.write("%s\n" % "\t".join(out_line))


if __name__ == "__main__":
  main()