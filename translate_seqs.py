#!/usr/bin/env python

'''
By Sean Lovett
8/2/2016
'''

from __future__ import print_function
from subprocess import Popen, PIPE
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys, optparse, os

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input ffn. [None, REQD]')

  opts, args = p.parse_args()

  seq_dict = translate_seqs(opts.input)
  
  write_fasta_from_dict(seq_dict, "%s.faa" % ".".join(opts.input.split(".")[:-1]))


def read_fasta_dict(file):
  names, seqs = read_fasta_lists(file)
  fasta_dict = dict(zip(names, seqs))
  return fasta_dict

def read_fasta_lists(file):
  fin = open(file, 'r')
  count=0

  names=[]
  seqs=[]
  seq=''
  for line in fin:
    line=line.strip()
    if line and line[0] == '>':                #indicates the name of the sequence
      count+=1
      names.append(line[1:].strip())
      if count>1:
        seqs.append(seq.strip())
      seq=''
    else:
      seq += line
  seqs.append(seq)

  return names, seqs

def write_fasta_from_dict(seq_dict, out_file):
  with open(out_file,"w") as fout:
    for key, seq in seq_dict.items():
      fout.write(">%s\n%s\n" % (key, seq))

def translate_seqs(in_fasta):
  seq_dict = read_fasta_dict(in_fasta)
  for key, seq in seq_dict.items():
    seq_dict[key] = translate_seq(seq)
  return seq_dict


# def translate_seq(in_seq):
#   code={'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C', 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 'TTA': 'L', 'TCA': 'S', 'TAA': 'STOP', 'TGA': 'STOP', 'TTG': 'L', 'TCG': 'S', 'TAG': 'STOP', 'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
#   
#   translated_list = []
#   for i in range(0,len(in_seq),3):
#     codon = in_seq[i:i+3]
#     if "N" in codon:
#       translated_list.append("X")
#     elif len(codon) < 3:
#       break
#     elif code[codon] == "STOP":
#       translated_list.append("*")
#     else:
#       translated_list.append(code[codon])
#   return "".join(translated_list)


def translate_seq(in_seq):
  if not len(in_seq) % 3 == 0:  
    for i in range(3-(len(in_seq)%3)):
      in_seq += "N"
  return str(Seq(in_seq, generic_dna).translate())


if __name__ == "__main__":
  main()