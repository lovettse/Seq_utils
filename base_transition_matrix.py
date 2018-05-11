#!/usr/bin/env python

from __future__ import division
from numpy import std
from decimal import Decimal
import optparse


def main():
  #To parse command line
  usage = "usage: %prog -i input.txt -r ref.fasta -g ref.gff [options] > commandlog.txt"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input pileup [None,REQD]')
  p.add_option('-o', '--output', help='Output [None,REQD]')
  
  opts, args = p.parse_args()

  bases = ['A','C','G','T']
  purines = ['A', 'G']
  pyrimidines = ['C', 'T']
  
  out_matrix = {}

  for base1 in bases:
    out_matrix[base1] = {}
    for base2 in bases:
      out_matrix[base1][base2] = [0,0]

  with open(opts.input, "r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      ref_base = cols[2].upper()
      pile_bases = cols[4].upper()
      prev_base = ""
      
      for pile_base in pile_bases:
        if prev_base.isdigit():
          continue
        # If the base is reference, increment numerator of ref_base/ref_base
        if pile_base in (".",","):
          out_matrix[ref_base][ref_base][0] += 1
        # If base is not reference, but, is a base (not read start/end, indel), increment numerator of base/ref_base 
        elif pile_base in bases:
          out_matrix[ref_base][pile_base][0] += 1
        if pile_base in bases or pile_base in pile_base in (".",","):
          for base in bases:
            # Increment denominator each base in ref_base row
            out_matrix[ref_base][base][1] += 1
        prev_base = pile_base

  transition_count = out_matrix['A']['G'][0] + out_matrix['G']['A'][0] + out_matrix['C']['T'][0] + out_matrix['T']['C'][0]
  transversion_count = out_matrix['A']['C'][0] + out_matrix['C']['A'][0] + out_matrix['G']['T'][0] + out_matrix['T']['G'][0] + out_matrix['A']['T'][0] + out_matrix['T']['A'][0] + out_matrix['G']['C'][0] + out_matrix['C']['G'][0]
  tstv = transition_count/transversion_count


  with open(opts.output,"w") as fout:
    fout.write("\t%s\n" % "\t".join(bases))
    for base1 in bases:
      out_line = list(base1)
      for base2 in bases:
        freq = (out_matrix[base1][base2][0])/(out_matrix[base1][base2][1])
        #out_line.append("%s/%s" % (out_matrix[base1][base2][0],out_matrix[base1][base2][1]))
        out_line.append("{:.4E}".format(Decimal(freq)))
      fout.write("%s\n" % "\t".join(out_line))


  freq_dict = {}
  avg_printed = 0
  for base1 in bases:
    for base2 in bases:
      if base1 != base2:
        freq_dict[(out_matrix[base1][base2][0])/(out_matrix[base1][base2][1])] = (base1, base2)
  avg_freq = sum(freq_dict)/len(freq_dict)
  with open(opts.input+".list", "w") as fout:
    fout.write("Mutation\tType\tFrequency\tFraction\n")
    fout.write("Total\t\t%s\t1\n" %("{:.4E}".format(Decimal(sum(freq_dict)))))
    for freq in sorted(freq_dict, reverse=True):
      if (freq_dict[freq][0] in purines and freq_dict[freq][1] in purines) or (freq_dict[freq][0] in pyrimidines and freq_dict[freq][1] in pyrimidines):
        mut_type = "trs"
      elif (freq_dict[freq][0] in purines and freq_dict[freq][1] in pyrimidines) or (freq_dict[freq][0] in pyrimidines and freq_dict[freq][1] in purines):
        mut_type = "trv"
      else:
        print "Invalid mut_type"
        exit()
      if freq < avg_freq and not avg_printed:
        fout.write("Average\t\t%s\t%0.4f\n" % ("{:.4E}".format(Decimal(avg_freq)), avg_freq/sum(freq_dict)))
        avg_printed += 1
      fout.write("%s -> %s\t%s\t%s\t%0.4f\n" % (freq_dict[freq][0], freq_dict[freq][1], mut_type, "{:.4E}".format(Decimal(freq)), freq/sum(freq_dict)))
#    fout.write("Transitions: %s\n" % transition_count)
#    fout.write("Transversions: %s\n" % transversion_count)
#    fout.write("ts/tv: %s\n" % tstv)
          

if __name__ == '__main__':
  main()