#!/usr/bin/env python

import optparse
from sets import Set

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input genedump file from single cell workflow[None, REQD]')
  p.add_option('-o', '--output', help='Output UMI expression matrix [None,REQD]')

  opts, args = p.parse_args()
  
  # First key = UMI, Second key = gene, Value = count
  umi_dict = {}
  # Key = gene, value = list of coexpression values for genes of interest
  out_dict = {}
  # Set of genes with any reads assigned to any UMI
  genes = Set()


  with open(opts.input,"r") as fin:
    fin.readline()
    for line in fin:
      cols = line.strip().split(",")
      barcode = cols[0]
      gene = cols[1]
      count = int(cols[2])
      
      if not barcode in umi_dict:
        umi_dict[barcode] = {}
      umi_dict[barcode][gene] = count
      
      genes.add(gene)

  with open(opts.output, "w") as fout:
    fout.write("%s\n" % ",".join(sorted(genes)))
    for umi in umi_dict:
      fout.write("%s" % umi)
      for gene in sorted(genes):
        if gene in umi_dict[umi]:
          fout.write(",%s" % (umi_dict[umi][gene]))
        else:
          fout.write(",0")
      fout.write("\n")


if __name__=="__main__":
  main()