#!/usr/bin/env python

'''
Sean Lovett 2017-10-3

Assigns frequencies to each potential base->base change.
Does not consider indels at all.
Does not use coverage information provided in pileup. Instead counts bases that are not indels.

Would love to have this handle indels as well, but that's not happening right now.
'''

from __future__ import division
from operator import itemgetter
from decimal import Decimal
from numpy import std
import optparse

bases = ('A','C','G','T')
purines = ('A', 'G')
pyrimidines = ('C', 'T')
indels = ('+', '-')


def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input pileup [None,REQD]')
  p.add_option('-o', '--output', help='Output [None,REQD]')
  p.add_option('-f', '--freqcut', type='float', default= 0.02, help = 'Threshold frequency. Frequency values below this will be assumed to be zero [0.02]')
  p.add_option('-d', '--depthcut', type='int', default = 200, help= 'Threshold depth. Bases with depth below this value will not be considered [200]')
  p.add_option('-r', '--regions', help= "File describing regions to analyze, should have 3 columns (tab-delimited) - contig name, start, stop [None]")
  p.add_option('--indels', default=False, action="store_true", help="Turn on indel calling. DON'T TRUST THIS YET. [False]")
  
  opts, args = p.parse_args()

  if opts.regions:
    regions = read_region_file(opts.regions)
  else:
    regions = 0

  if not opts.indels:
    global indels
    indels = ()

  # out_matrix is a dict of dicts
  # First key is the reference base, second key is the aligned base
  # At first, out_matrix will be the sum of frequencies at which a given base is a aligned to a given reference base
  # When all bases have been counted, those values will be divided by the base counts to give an average frequency
  out_matrix = get_freq_matrix(opts.input, bases, opts.freqcut, opts.depthcut, regions)

  # Prints out_matrix to file
  print_out_matrix(opts.output, out_matrix)

  # Print list output to file
  print_list_output(opts.output,out_matrix)


def read_region_file(region_file):
  regions = {}

  with open(region_file) as fin:
    for line in fin:
      cols = line.strip().split("\t")
      if not cols[0] in regions:
        regions[cols[0]] = []
      regions[cols[0]].append((int(cols[1]), int(cols[2])))

  return regions


def print_list_output(out_list, out_matrix):
  freq_list = []
  avg_printed = 0
  for ref_base in bases:
    for aligned_base in bases+indels:
      if ref_base != aligned_base:
        freq_list.append((ref_base, aligned_base, out_matrix[ref_base][aligned_base]))
  freq_list.sort(key=itemgetter(2), reverse=True)
  
  avg_freq = sum([freq[2] for freq in freq_list])/len(freq_list)
  with open(out_list+".list", "w") as fout:
    fout.write("Mutation\tType\tFrequency\tFraction\n")
    fout.write("Total\t\t%0.4f\t1\n" % (sum([frequ[2] for frequ in freq_list])))
    for freq in freq_list:
      if (freq[0] in purines and freq[1] in purines) or (freq[0] in pyrimidines and freq[1] in pyrimidines):
        mut_type = "trs"
      elif (freq[0] in purines and freq[1] in pyrimidines) or (freq[0] in pyrimidines and freq[1] in purines):
        mut_type = "trv"
      else:
        mut_type= "indel"
      if freq[2] < avg_freq and not avg_printed:
        fout.write("Average\t\t%s\t%0.4f\n" % (avg_freq, avg_freq/sum([frequ[2] for frequ in freq_list])))
        avg_printed += 1
      if freq[2] > 0:
        fout.write("%s -> %s\t%s\t%s\t%0.4f\n" % (freq[0], freq[1], mut_type, freq[2], freq[2]/sum([frequ[2] for frequ in freq_list])))
      else:
        fout.write("%s -> %s\t%s\t0\t0\n" % (freq[0], freq[1], mut_type))


def print_out_matrix(out_mat_file, out_matrix):
  with open(out_mat_file,"w") as fout:
    fout.write("\t%s\n" % "\t".join(bases+indels))
    for ref_base in bases:
      out_line = list(ref_base)
      for aligned_base in bases+indels:
        #out_line.append(str("{:.4E}".format(Decimal(out_matrix[ref_base][aligned_base]))))
        out_line.append( "%0.4f" % (out_matrix[ref_base][aligned_base]))
      fout.write("%s\n" % "\t".join(out_line))


def get_freq_matrix(in_pile, bases, freqcut, depthcut, regions):
  out_matrix = {}
  
  base_counts = {}
  for base in bases:
    base_counts[base] = 0

  for ref_base in bases:
    out_matrix[ref_base] = {}
    for aligned_base in bases+indels:
      out_matrix[ref_base][aligned_base] = 0

  with open(in_pile, "r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      ref_base = cols[2].upper()
      pile_bases = iter(cols[4].upper())
      prev_base = ""
      
      contig_name = cols[0]
      position = int(cols[1])

      base_depth = 0

      # Check that the position is included in the list of user-provided regions
      valid_position = 0
      # If regions were not provided, all positions are valid
      if not regions:
        valid_position = 1
      # If regions were provided and the contig_name appears at least once, scan the list of regions for the contig
      elif regions and contig_name in regions:
        for region in regions[contig_name]:
          if position >= region[0] and position <= region[1]:
            valid_position = 1
      
      if valid_position:
        counts = {}
        for base in bases+indels:
          counts[base] = 0
        
        for pile_base in pile_bases:
          if prev_base.isdigit():
            continue
          # If the base is reference, add it count for ref_base
          if pile_base in (".",","):
            base_depth += 1
            counts[ref_base] += 1
          # If base is not reference, but is a base (not read start/end, indel), increment count for that base
          elif pile_base in bases+indels:
            base_depth += 1
            counts[pile_base] += 1
            if pile_base in indels:
              print pile_base + pile_bases.next() + pile_bases.next()
          # elif pile_base == "+":
          #   pass
          # elif pile_base == "-":
          #   pass
          
          prev_base = pile_base
  
        if base_depth >= depthcut:
          base_counts[ref_base] += 1
          # Add frequency to the output_matrix
          for aligned_base in bases+indels:
            if counts[aligned_base]/base_depth > freqcut or aligned_base in indels:
              out_matrix[ref_base][aligned_base] += counts[aligned_base]/base_depth

  # Go through matrix and divide the frequencies by the number of occurences of the reference base to get average frequency
  for ref_base in bases:
    for aligned_base in bases+indels:
      frequency = out_matrix[ref_base][aligned_base]
      # In the event that there are none of a given base in the provided region, just set frequency to zero
      if base_counts[ref_base] == 0:
        avg_frequency = 0
      else:
        avg_frequency = frequency/base_counts[ref_base]
      out_matrix[ref_base][aligned_base] = avg_frequency
      
  return out_matrix


if __name__ == '__main__':
  main()
