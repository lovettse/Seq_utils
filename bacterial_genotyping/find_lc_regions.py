#!/usr/bin/env python

'''
By Sean Lovett
8/19/2016

This script identifies low-complexity regions using the "dust" method in prinseq.

Its primary output is a bed file describing regions that have been filtered.

It depends on the cut_ref_fastq.py script and prinseq. 
'''

import optparse
from subprocess import Popen, PIPE

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  #Input/output files
  p.add_option('-i', '--input', help='Input fastq created by cut_ref_fastq.py [None,REQD]')
  p.add_option('-p', '--prinseq', default="/data/_pipeline_software/prinseq-lite-0.20.3/prinseq-lite.pl", help='How to call prinseq [/data/_pipeline_software/prinseq-lite-0.20.3/prinseq-lite.pl]')
  p.add_option('-f', '--filtered', help="Optional starting place. Fastq that has already been filtered with prinseq [None]")
  p.add_option('-o', '--output', help='Output bed. [`basename opts.input .fastq`_lc_regions.bed]')
  p.add_option('-t', '--threshold', default=3, type="int", help='lc_threshold parameter for prinseq [3]')
  
  opts, args = p.parse_args()
  
  if not opts.output:
    opts.output = "%s_lc_regions.bed" % ".".join(opts.input.split(".")[:-1])
  
  headers = []
  filtered_bases = {}

  if opts.filtered:
    filtered_fastq = opts.filtered
  else:
    filtered_fastq = run_prinseq(opts)

  with open(filtered_fastq,"r") as fin:
    for line in fin:
      if line.startswith("@"):
        headers.append(line.strip()[1:])

  for header in headers:
    chrom = "_".join(header.split("_")[:-1])
    pos = int(header.split("_")[-1])
    if chrom in filtered_bases:
      filtered_bases[chrom].append(pos)
    else:
      filtered_bases[chrom] = []
      filtered_bases[chrom].append(pos)

  make_out_list(filtered_bases,opts)
  make_out_bed(filtered_bases,opts)


def make_out_list(filtered_bases,opts):
  out_name = "%s_lc_bases.txt" % ".".join(opts.output.split(".")[:-1])

  with open(out_name,"w") as fout:
    for chrom, pos_list in filtered_bases.iteritems():
      for pos in sorted(pos_list):
        fout.write("%s\t%s\n" % (chrom, pos))


def make_out_bed(filtered_bases, opts):
  with open(opts.output,"w") as fout:
    for chrom, pos_list in filtered_bases.iteritems():
      start = 0
      stop = 0
      for pos in sorted(pos_list):
        if not start and not stop:
          start = pos
          stop = pos
        elif pos-stop == 1:
          stop=pos
        else:
          fout.write("%s\t%d\t%d\n" % (chrom, start-1, stop))
          start = 0
          stop = 0

  
def run_prinseq(opts):
  out_name = "%s_lc_filtered" % ("".join(opts.input.split(".")[:-1]))
  cline = "%s -fastq %s -out_bad %s -out_good null -lc_method dust -lc_threshold %s 2> prinseq.err" % (opts.prinseq, opts.input, out_name, opts.threshold)
  print(cline)
  prinseq_run=Popen(cline, shell=True, stdout=PIPE, stderr=PIPE)
  prinseq_run.wait()
  return "%s.fastq" % out_name


if __name__ == "__main__":
  main()
