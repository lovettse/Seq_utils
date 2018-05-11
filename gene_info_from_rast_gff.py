#!/usr/bin/env python

from __future__ import division
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input gff.\nThis definitely works with gffs made by RAST. It might work with others [None, REQD]')
  p.add_option('-o', '--output', help="Output gene info [None, REQD]")

  opts, args = p.parse_args()

  info_dict = {}
  
  with open(opts.input, "r") as fin:
    for line in fin:
      if line.strip() == "##FASTA":
        break
      if not line.startswith("#"):
        cols = line.strip().split("\t")
        if cols[2] == "CDS":
          gene_id = cols[8].split("=")[1].split(";")[0]
          chrom = cols[0]
          if "Name" in cols[8].split(";")[1]:
            gene_name = cols[8].split(";")[1].split("=")[1]
          else:
            gene_name = ""
          start = cols[3]
          end = cols[4]
          info_dict[gene_id] = (chrom,gene_name,start,end)

  with open(opts.output, "w") as fout:
    fout.write('Feature_ID\tChrom_ID\tGene_Name\tStart\tEnd\n')
    for gene_id in sorted(info_dict):
      fout.write("%s\t%s\n" % (gene_id, "\t".join(info_dict[gene_id])))

  
if __name__=="__main__":
  main()