#!/usr/bin/env python

from __future__ import division
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input htseq file[None, REQD]')
  p.add_option('-o', '--output', default="/dev/stdout", help='Output summary [stdout]')

  opts, args = p.parse_args()

  htseq_dict = {}
  out_dict = {}
  total_reads = 0
  gene_reads = 0
  
  with open(opts.input, "r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      htseq_dict[cols[0]] = int(cols[1])

  for gene in htseq_dict:
    total_reads += htseq_dict[gene]
    if not gene.startswith("__"):
      gene_reads += htseq_dict[gene]
  
  for gene in htseq_dict:
    if gene.startswith("__"):
      category = gene[2:]
      out_dict[category] = htseq_dict[gene]/total_reads

  with open(opts.output, "w") as fout:
    fout.write("total_reads\t%s\n" % (total_reads))
    fout.write("gene_read_count\t%s\n" % (gene_reads))
    fout.write("gene_read_perc\t%0.2f%%\n" % (gene_reads/total_reads*100))
    for line in sorted(out_dict):
      fout.write("%s\t%0.2f%%\n" % (line, out_dict[line]*100))    

if __name__=="__main__":
  main()