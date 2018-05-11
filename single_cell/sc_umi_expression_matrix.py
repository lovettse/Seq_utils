#!/usr/bin/env python

import optparse
from sets import Set

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-i', '--input', help='Input genedump file from single cell workflow[None, REQD]')
  p.add_option('-l', '--list', help='Input list of genes of interest, will use all genes if not given [None]')
  p.add_option('-o', '--output', help='Output UMI expression matrix [None,REQD]')
  p.add_option('-c', '--count', dest="count", action="store_true", help='Output read count matrix, default is to count UMIs with any level of co-expression equally[False]' )

  opts, args = p.parse_args()
  
  # First key = UMI, Second key = gene, Value = count
  umi_dict = {}
  # Key = gene, value = list of coexpression values for genes of interest
  out_dict = {}
  # Set of genes with any reads assigned to any UMI
  genes = Set()
  # User-provided set of genes of interest
  genes_of_interest = []

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

  if opts.list:
    with open(opts.list, "r") as fin:
      for line in fin:
        genes_of_interest.append(line.strip())
    for gene in genes_of_interest:
      if not gene in genes:
        print "%s was not found in the results file. Either it was not expressed or you are using an incorrect gene name" % gene
  else:
    genes_of_interest = genes

  for gene in genes:
    out_dict[gene] = []
    total = 0
    for gene_of_interest in genes_of_interest:
      co_expressed = 0
      for barcode in umi_dict:
        if gene_of_interest in umi_dict[barcode] and gene in umi_dict[barcode]:
          if opts.count:
            co_expressed += umi_dict[barcode][gene]
          else:
            co_expressed += 1
      total += co_expressed
      out_dict[gene].append(str(co_expressed))
    out_dict[gene].append(str(total))

  with open(opts.output, "w") as fout:
    fout.write("gene\t%s\ttotal\n" % "\t".join(genes_of_interest))
    for gene in sorted(out_dict):
      fout.write("%s\t%s\n" % (gene, "\t".join(out_dict[gene])))

if __name__=="__main__":
  main()