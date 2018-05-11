#!/usr/bin/python

from Bio import SeqIO
import optparse, os

usage = '%prog [options]'
p = optparse.OptionParser()

p.add_option('-i', '--input', help='Input fasta [None, REQD]')
p.add_option('-o', '--output', help='Output fasta [None, REQD]')
p.add_option('-a', '--all', action="store_true", default=False, help='Include everything with a gene annotation, not just CDS')


opts, args = p.parse_args()

gbk_filename = opts.input
faa_filename = opts.output


input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    output_handle.write("Feature_ID\tChrom_ID\tGene_Name\tStart\tEnd\n")
    if not opts.all:
        for seq_feature in seq_record.features:
            start = int(str(seq_feature.location.start).replace(">","").replace("<","")) + 1
            end = int(str(seq_feature.location.end).replace(">","").replace("<",""))
            if seq_feature.type=="CDS":
                if "gene" in seq_feature.qualifiers:
                    gene = seq_feature.qualifiers['gene'][0]
                else:
                    gene = ""
                output_handle.write("%s\t%s\t%s\t%s\t%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_record.id,
                       gene,
                       start,
                       end))
    else:
        for seq_feature in seq_record.features:
            start = int(str(seq_feature.location.start).replace(">","").replace("<","")) + 1
            end = int(str(seq_feature.location.end).replace(">","").replace("<",""))
            if seq_feature.type=="gene":
                if "gene" in seq_feature.qualifiers:
                    gene = seq_feature.qualifiers['gene'][0]
                else:
                    gene = ""
                output_handle.write("%s\t%s\t%s\t%s\t%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_record.id,
                       gene,
                       start,
                       end))

output_handle.close()
input_handle.close()
