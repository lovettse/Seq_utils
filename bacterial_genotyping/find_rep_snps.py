#!/usr/bin/env python

'''
By Sean Lovett
8/18/2016
'''

from __future__ import print_function
import optparse

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input bed file created by rep_regions_bybase.py [None,REQD]')
  p.add_option('-l', '--list', help='Input list of snps, one per line, tab delimited. Chrom\tpos [None,REQD]')
  p.add_option('-o', '--output', help='Output file. Two columns, tab delimited. Position\t(bool)repetitive [None, REQD]')
  p.add_option('-w', '--window', type="int", default=50, help='Window size used in rep_regions_bybase.py [50]')
  # p.add_option('-p', '--pad', type="int", default=50, help='"Pad" distance around repetitive regions [50]')

  opts, args = p.parse_args()
  
  snp_pos_list = []
  rep_region_list = []
  
  with open(opts.list,"r") as fin:
    for line in fin:
      if not line.startswith("#") and not line.startswith("CHROM"):
        cols = line.strip().split("\t")
        chrom = cols[0]
        pos = int(cols[1])
        snp_pos_list.append((chrom,pos))

  with open(opts.input,"r") as fin:
    for line in fin:
      cols = line.strip().split("\t")
      chrom = cols[0]
      begin = int(cols[1])
      end = int(cols[2]) + opts.window
      
      rep_region_list.append((chrom,begin,end))

  with open(opts.output,"w") as fout:
    for snp_pos in snp_pos_list:
      rep = False
      for rep_region in rep_region_list:
        if snp_pos[0] == rep_region[0] and snp_pos[1] >= rep_region[1] and snp_pos[1] <= rep_region[2]:
          fout.write("%s\t%s\tTrue\n" % (snp_pos[0],snp_pos[1]))
          rep = True
          break
      # if not rep:
      #   for rep_region in rep_region_list:
      #     if snp_pos[0] == rep_region[0] and snp_pos[1] >= rep_region[1] - opts.pad and snp_pos[1] <= rep_region[2] + opts.pad:
      #       fout.write("%s\t%s\tWithin %s\n" % (snp_pos[0],snp_pos[1],opts.pad))
      #       rep = True
      #       break
      #   if not rep:
      #     fout.write("%s\t%s\tFalse\n" % (snp_pos[0],snp_pos[1]))
      if not rep:
        # I hate this magic number. Show me a genome with greater than 10MB between any two repeats and I'll change it.
        min_dist = 10000000
        for rep_region in rep_region_list:
          if rep_region[0] == snp_pos[0]:
            if abs(rep_region[1]-snp_pos[1]) < min_dist:
              min_dist = abs(rep_region[1]-snp_pos[1])
            if abs(rep_region[2]-snp_pos[1]) < min_dist:
              min_dist = abs(rep_region[2]-snp_pos[1])
        fout.write("%s\t%s\t%s\n" % (snp_pos[0], snp_pos[1], min_dist))


if __name__ == "__main__":
  main()