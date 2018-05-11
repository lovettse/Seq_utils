#!/usr/bin/env python

import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input genedump.csv [None,REQD]')
  p.add_option('-m', '--meta', help='Input metadata [None,REQD]')
  p.add_option('-g', '--gtf', help='Input GTF used to make STAR index [None,REQD]')
  p.add_option('-t', '--term', default="gene", help='Term to look for in GTF. Should be column 2 of lines that have gene_name. Change to CDS for non-Ensembl GTFs[gene]')
  p.add_option('-o', '--output', help='Output genedump with cell type info [None,REQD]')
  
  opts, args = p.parse_args()
  
  meta_dict = {}
  genename_dict = {}
  
  with open(opts.gtf, "r") as fin:
    for line in fin:
      if not line.startswith("#"):
        cols = line.strip().split("\t")
        if cols[2] == opts.term:
          info = cols[8].split(";")
          geneid = info[0].split(" ")[1][1:-1]
          genename = info[2].split(" ")[2][1:-1]
          genename_dict[geneid] = genename
    

  with open(opts.meta, "r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      meta_dict[cols[0]] = cols[2]

  with open(opts.input, "r") as fin:
    with open(opts.output,"w") as fout:
      fin.readline()
      fout.write("Barcode,Cell_type,Gene,Gene_name,Count\n")
      for line in fin:
        cols = line.strip().split(",")
        if cols[1] in genename_dict:
          genename = genename_dict[cols[1]]
        else:
          genename= "N/A"
        fout.write("%s,%s,%s,%s,%s\n" % (cols[0], meta_dict[cols[0]], cols[1], genename,cols[2]))
      

if __name__ == '__main__':
  main()