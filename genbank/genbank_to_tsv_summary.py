#!/usr/bin/python

from Bio import SeqIO
import optparse, os

usage = '%prog [options]'
p = optparse.OptionParser()

p.add_option('-i', '--input', help='Input Genbank file [None, REQD]')
p.add_option('-o', '--output', help='Output text file [None, REQD]')

opts, args = p.parse_args()

gbk_filename = opts.input
faa_filename = opts.output

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

output_handle.write("locus_tag\tgene\tcontig\tstart\tend\tstrand\ttype\tproduct\tnuc_seq\taa_seq\n")

for seq_record in SeqIO.parse(input_handle, "genbank") :
  print "Dealing with GenBank record %s" % seq_record.id
  last_was_gene = 0
  first_gene = 0
  for seq_feature in seq_record.features :
    if seq_feature.type=="gene":
      first_gene = 1
      if 'gene' in seq_feature.qualifiers:
        gene = seq_feature.qualifiers['gene'][0]
      else:
        gene = ""
      if seq_feature.strand < 0:
        strand = "-"
      elif seq_feature.strand > 0:
        strand = "+"
      output_handle.write("%s\t%s\t%s\t%s\t%s\t%s" % (seq_feature.qualifiers['locus_tag'][0],
                                                  gene,
                                                  seq_record.id,
                                                  seq_feature.location.start+1,
                                                  seq_feature.location.end+1,
                                                  strand))
      last_was_gene = 1
    elif seq_feature.type=="CDS" and last_was_gene:
      if "translation" in seq_feature.qualifiers:
        trans = seq_feature.qualifiers['translation'][0]
      else:
        trans = ""
      if "pseudo" in seq_feature.qualifiers:
        output_handle.write("\tCDS (pseudo)\t%s\t%s\t%s\n" % (seq_feature.qualifiers['product'][0],
                                                              seq_feature.extract(seq_record.seq),
                                                              trans))
      else:
        output_handle.write("\tCDS\t%s\t%s\t%s\n" % (seq_feature.qualifiers['product'][0],
                                                     seq_feature.extract(seq_record.seq),
                                                     trans))
      last_was_gene = 0
    elif first_gene and last_was_gene:
      if 'product' in seq_feature.qualifiers:
        product = seq_feature.qualifiers['product'][0]
      else:
        product = ''
      output_handle.write("\t%s\t%s\t%s\n" % (seq_feature.type,
                                              product,
                                              seq_feature.extract(seq_record.seq)))
      last_was_gene = 0

output_handle.close()
input_handle.close()