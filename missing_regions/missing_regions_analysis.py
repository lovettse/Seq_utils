#!/usr/bin/env python

from __future__ import division
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)

  #Input/output files
  p.add_option('-s', '--samples', help='Input list of all samples. [None, REQD]')
  p.add_option('-c', '--clade', help='Input list of samples within clade of interest [None, REQD]')
  p.add_option('-g', '--genes', help='Input list of all genes present in reference [None,REQD]')
  p.add_option('-f', '--gff', help='Input gff annotating reference gene set [None,REQD]')
  p.add_option('-o', '--output', help='Output list of genes specific to clade [None,REQD]')
  p.add_option('-t', '--table', help='Output summary table [None,REQD]')

  opts, args = p.parse_args()

  samples = set()
  clade = set()
  genes = set()
  missing_genes = {}
  
  missing_from_nonclade = {}
  
  with open(opts.samples,"r") as fin:
    for line in fin:
      samples.add(line.strip())
  with open(opts.clade,"r") as fin:
    for line in fin:
      clade.add(line.strip())
  with open(opts.genes,"r") as fin:
    for line in fin:
      genes.add(line.strip())

  core_genome = set(genes)
  clade_genome = set(genes)
  nonclade_genome = set(genes)
  nonclade = set(samples).difference(set(clade))
  missing_from_nonclade = set()
  present_in_clade = set(genes)

  for sample in samples:
    missing_genes[sample] = set()
    with open("%s_missing_genes.txt" % sample, "r") as fin:
      fin.readline()
      for line in fin:
        cols = line.strip().split("\t")
        missing_genes[sample].add(cols[0])

  for gene in genes:
    missing_count = 0
    for sample in nonclade:
      if gene in missing_genes[sample]:
        nonclade_genome.discard(gene)
        missing_count += 1
        core_genome.discard(gene)
    if missing_count >= len(nonclade):
      missing_from_nonclade.add(gene)
    for sample in clade:
      if gene in missing_genes[sample]:
        present_in_clade.discard(gene)
        core_genome.discard(gene)

  clade_specific_genome = missing_from_nonclade.intersection(present_in_clade)

  print "reference genome gene count: %s" % len(genes)
  print "core genome gene count: %s" % len(core_genome)
  print "clade core genome gene count: %s" % len(clade_genome)
  #print "nonclade genome gene count: %s" % len(nonclade_genome)
  print "clade specific gene count: %s" % len(clade_specific_genome)
  
  with open(opts.output,"w") as fout:
    with open(opts.gff, "r") as fin:
      for line in fin:
        if line.strip() == "##FASTA":
          break
        if not line.startswith("#"):
          cols = line.strip().split("\t")
          if cols[2] == "CDS":
            for info in cols[8].split(";"):
              if info.startswith("ID="):
                gene_id = info.split("=")[1]
              elif info.startswith("Name="):
                product = info.split("=")[1]
            if gene_id in clade_specific_genome:
              fout.write("%s\t%s\n" % (gene_id,product))

  if opts.table:
    with open(opts.table,"w") as fout:
      fout.write("gene\t%s\n" % "\t".join(samples))
      for gene in genes:
        fout.write("%s" % gene)
        for sample in samples:
          if gene in missing_genes[sample]:
            fout.write("\t0")
          else:
            fout.write("\t1")
        fout.write("\n")

if __name__=="__main__":
  main()