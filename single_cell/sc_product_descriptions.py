#!/usr/bin/env python

import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input file from single cell workflow[None, REQD]')
  p.add_option('-d', '--desc', help='Description file [None,REQD')
  p.add_option('-o', '--output', help='Output file [None,REQD]')

  opts, args = p.parse_args()
  
  count_dict = {}
  desc_dict = {}
  
  with open(opts.input,"r") as fin:
    for line in fin:
      if not line.startswith("#") and not line.startswith("\t"):
        cols = line.strip().split("\t")
        gene = cols[0]
        count = cols[6]
        count_dict[gene] = count

  with open(opts.desc,"r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      gene = cols[0]
      desc = cols[1]
      desc_dict[gene] = desc

  out_lines = []
  for gene in sorted(count_dict):
    if gene in desc_dict:
      out_lines.append("%s\t%s\t%s\n" % (gene, desc_dict[gene], count_dict[gene]))
    else:
      out_lines.append("%s\t\t%s\n" % (gene, count_dict[gene]))

  with open(opts.output, "w") as fout:
    for line in out_lines:
      fout.write(line)
    
    

if __name__=="__main__":
  main()