#!/usr/bin/python

from Bio import SeqIO
import optparse, os

usage = '%prog [options]'
p = optparse.OptionParser()

p.add_option('-i', '--input', help='Input fasta [None, REQD]')
p.add_option('-o', '--output', help='Output fasta [None, REQD]')

opts, args = p.parse_args()

gbk_filename = opts.input
faa_filename = opts.output


input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print "Dealing with GenBank record %s" % seq_record.id
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
           
            if 'translation' in seq_feature.qualifiers: 
                output_handle.write(">%s\n%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
