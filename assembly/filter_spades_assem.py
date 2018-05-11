#!/usr/bin/python

from __future__ import division
import optparse
from Bio import SeqIO

def main():


  usage = '%prog [options]'
  p = optparse.OptionParser()

  p.add_option('-i', '--input', help='Input assembly from SPAdes [None, REQD]')
  p.add_option('-l', '--length', default=200, type='int',help='Minimum contig length [200]')
  p.add_option('-c', '--cov', default=2, type='int', help='Minimum coverage [2]')
  p.add_option('-o', '--output', help="Output fasta [None, REQD]")

  opts, args = p.parse_args()
  
  out_seqs = []
  
  with open(opts.input, "r") as fin:
    for seq_record in SeqIO.parse(fin, "fasta"):
      info = seq_record.id.split("_")
      if int(info[3]) > opts.length and float(info[5]) > opts.cov:
        out_seqs.append((seq_record.id,seq_record.seq))
  
  with open(opts.output, "w") as fout:
    for out_seq in out_seqs:
      fout.write(">%s\n%s\n" % (out_seq[0], out_seq[1]))
    
  
if __name__=="__main__":
  main()