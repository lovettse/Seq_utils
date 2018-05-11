#!/usr/bin/env python

from Bio import SeqIO
from collections import Counter
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input viral multiple alignment [None,REQD]')
  p.add_option('-n', '--num', default=10, type= 'int', help='Number of bases required to agree with reference at beginning and end [10]')
  p.add_option('-o', '--output', help='Output trimmed alignment (not in same order as input) [None,REQD]')
  
  opts, args = p.parse_args()

  consensus = get_consensus(opts.input)

  write_fasta({opts.input+"_consensus": consensus}, "consensus.fasta")

  out_seqs = trim_seqs(opts.input, consensus, opts.num)
  
  write_fasta(out_seqs, opts.output)
  
  get_consensus(opts.output)
  


def trim_seqs(in_fasta, consensus, num):
  out_seqs = {}

  with open(in_fasta, "r") as fin:
    for record in SeqIO.parse(fin,"fasta"):
      header = record.description
      seq = record.seq.upper()
      rev_seq = seq[::-1]

      trim_seq_forw,count_forw = trim_seq(seq, consensus, num)
      trim_seq_rev,count_rev = trim_seq(rev_seq, consensus[::-1], num)

      out_seqs[record.description] = combine_seqs(trim_seq_forw, trim_seq_rev[::-1])
      # if count_forw > 30:
      #   print "%s\t%s" % (header, count_forw)
  return out_seqs


def combine_seqs(forw, rev):
  combined = ""
  counter = 0
  if len(forw) != len(rev):
    print "Seqs are different lengths!!!!!"
    exit
  for i in range(len(forw)):
    if forw[i] == "N" or rev[i] == "N":
      combined += "N"
      if forw[i] != "N" or rev[i] != "N":
        counter += 1
    elif forw[i] == rev[i]:
      combined += forw[i]
    else:
      print "Neither forward nor reverse is N, but they do not match"
      print "Forward: %s\nReverse: %s" % (forw[i], rev[i])
#  if counter:
#    print "%s Ns inserted" % (counter)
  return combined


def trim_seq(seq, consensus, base_count):
  new_seq = ""
  from_begin = 0
  first_notn = 0
  last_addedn = 0

  for i in range(len(seq)):
    if from_begin >= base_count:
      if seq[i] in ['A','C','G','T','-']:
        new_seq += seq[i]
      else:
        new_seq += "N"
    elif seq[i] not in ['A','C','G','T','-']:
      new_seq = ""
      for k in range(i+1):
        new_seq += "N"
      last_addedn = i
    else:
      if not first_notn: first_notn = i
      if seq[i] == consensus[i]:
        new_seq += seq[i]
        from_begin += 1
      else:
        new_seq = ""
        for j in range(i+1):
          new_seq += "N"
        last_addedn = i

  #print "%s Ns added" % (last_addedn - first_notn + 1)

  return new_seq, last_addedn - first_notn + 1


def write_fasta(seq_dict, out_fasta):
  with open(out_fasta, "w") as fout:
    for seq in seq_dict:
      fout.write(">%s\n%s\n" % (seq, seq_dict[seq]))


def get_consensus(in_fasta):
  consensus = ""
  seq_list = []
  snps = []

  with open(in_fasta, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      i = 0
      for base in record.seq.upper():
        if len(seq_list) == i:
          seq_list.append(Counter())
        if base in ['A','C','G','T','-']:
          seq_list[i].update(base)
        i+=1

  j = 0
  for base in seq_list:
    if base:
      consensus += base.most_common(1)[0][0]
    else:
      consensus += "N"
    del base['-']
    if len(base.most_common(5)) > 1:
      snps.append(base)
      print "%s\t%s" % (j+1, base)
    j += 1

  return consensus


if __name__ == '__main__':
  main()