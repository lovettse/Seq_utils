#!/usr/bin/env python

from Bio import SeqIO
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input _parsed.txt [None,REQD]')
  p.add_option('-f', '--fasta', help='Input fasta [None,REQD]')
  p.add_option('-o', '--output', help='Output fasta [None,REQD]')
  p.add_option('--query', default=False, action="store_true", help='Set if you\'d rather use the query name to extract [False]')
  
  opts, args = p.parse_args()

  names_to_extract = []
  out_fasta_dict = {}

  with open(opts.input,"r") as fin:
    header = fin.readline()
    for line in fin:
      cols = line.strip().split("\t")
      if not opts.query:
        subject_name = " ".join(cols[2].split(" ")[1:])
        names_to_extract.append(subject_name)
      elif opts.query:
        names_to_extract.append(cols[0])
      
      
  with open(opts.fasta, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      if record.description in names_to_extract:
        out_fasta_dict[record.description] = str(record.seq)

  with open(opts.output, "w") as fout:
    for seq in out_fasta_dict:
      fout.write(">%s\n%s\n" % (seq, out_fasta_dict[seq]))


if __name__ == '__main__':
  main()