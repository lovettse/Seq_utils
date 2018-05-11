#!/usr/bin/python

from Bio import SeqIO
import optparse, os

usage = '%prog [options]'
p = optparse.OptionParser()

p.add_option('-i', '--input', help='Input fasta [None, REQD]')
p.add_option('-o', '--output', help='Output fasta [None, REQD]')

opts, args = p.parse_args()

gbk_filename = opts.input
fna_filename = opts.output

input_handle  = open(gbk_filename, "r")
output_handle = open(fna_filename, "w")

print "Dealing with GenBank file %s" % gbk_filename
SeqIO.convert(gbk_filename, "genbank", fna_filename, "fasta")

output_handle.close()
input_handle.close()
