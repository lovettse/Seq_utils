#!/usr/bin/env python

from Bio import SeqIO
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fasta [None,REQD]')
  p.add_option('-r', '--remove', help="List of fasta headers to remove, must exactly match [None,REQD]")
  p.add_option('-o', '--output', help='Output fasta [None,REQD]')
  
  opts, args = p.parse_args()

  remove_list = []
  out_seqs = {}

  with open(opts.remove, "r") as fin:
    for line in fin:
      remove_list.append(line.strip())

  with open(opts.input, "r") as fin:
    for record in SeqIO.parse(fin,"fasta"):
      if not record.id in remove_list:
        out_seqs[record.id] = str(record.seq)

  with open(opts.output, "w") as fout:
    for seq in out_seqs:
      fout.write(">%s\n%s\n" % (seq, out_seqs[seq]))

if __name__ == '__main__':
  main()