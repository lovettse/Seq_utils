#!/usr/bin/env python

'''
By Sean Lovett
9/9/2016
'''

from __future__ import print_function
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input list of ffns, one per line. [None, REQD]')
  p.add_option('-o', '--output', default="existence_matrix.txt", help='Output matrix, tab-delimited. [existence_matrix.txt]')
  p.add_option('-l', '--list', help="List of gene names [None,REQD]")
  
  opts, args = p.parse_args()
  
  in_fastas = []
  existence_dict = {}
  gene_list = []
  output = []
  ortho_dict = {}
  
  with open(opts.input,"r") as fin:
    for line in fin:
      in_fastas.append(line.strip())
 
  for fasta in in_fastas:
    species = fasta.split("_")[0]
    existence_dict[species] = []
    with open(fasta,"r") as fin:
      for line in fin:
        if line.startswith(">"):
          existence_dict[species].append(line.strip())

  with open(opts.list,"r") as fin:
    for line in fin:
#      gene_list.append(line.strip().split("\t")[0].upper())
      gene_list.append(line.strip().upper())

  header = []
  header.append("gene")
  
  for species in existence_dict:
    header.append(species)

  header.append("count")
  
  for gene in sorted(gene_list):
    ortho_dict[gene] = []
    out_line = []
    out_line.append(gene)
    exists = 0
    for species in existence_dict:
      found = False
      for def_line in existence_dict[species]:
        if gene in def_line.upper().split("|"):
          ortho_dict[gene].append(def_line[1:])
          found = True
          exists += 1
          break
      if found:
        out_line.append("1")
      else:
        ortho_dict[gene].append("")
        out_line.append("0")
    out_line.append(str(exists))
    output.append(out_line)

  with open(opts.output, "w") as fout:
    fout.write("%s\n" % "\t".join(header))
    for line in output:
      fout.write("%s\n" % "\t".join(line))

  with open("ortho_table.txt", "w") as fout:
    for gene in sorted(ortho_dict):
      fout.write("%s\t%s\n" % (gene, "\t".join(ortho_dict[gene])))
    
  
if __name__ == "__main__":
  main()
