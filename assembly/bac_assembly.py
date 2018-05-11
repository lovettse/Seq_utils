#!/usr/bin/env python

'''
bac_assembly.py
Sean Lovett
11/10/2017

This script is intended to use PacBio or Nanopore reads to recreate the linear piece of DNA inserted into a circular BAC vector.

Inputs:
Subread fastq from PacBio (need to get some Nanopore data to try this out on and know which file to use)
Fasta containing the expected linear represention of the cloning vector (this may require you to rotate it around the site at which it was cut)

Outputs:
Fasta containing the sequence inserted into the BAC
Intermediate files including BLAST and Canu output

Dependencies:
Python 2 (I think either 2.6 or 2.7 will work)
Biopython
Java 1.8
Canu
BLAST
split_blast_multi_faster.py (or something that produces equivalent and equivalently named output)
trim_shift_BAC.py (internal script found in /data/_scripts/slovett on Prod1)

The paths for all external programs called can be set at runtime. The default values will work for Prod1.
'''

from __future__ import division
from Bio import SeqIO
from subprocess import Popen, PIPE
import optparse,sys

def main():
  #To parse command line
  usage = "usage: %prog [options]"
  p = optparse.OptionParser(usage)
  
  p.add_option('-i', '--input', help='Input fastq [None, REQD]')
  p.add_option('-n', '--name', default= "bac_assem", help='Assembly name [bac_assem]')
  p.add_option('-l', '--log', default = "command.log", help="Command log [command.log]")
  p.add_option('-b', '--backbone', help='Input backbone fasta [None,REQD]')
  p.add_option('-o', '--output', default= "trimmed_BAC.fasta", help='Output fasta [trimmed_BAC.fasta]')
  p.add_option('--nanopore', default=False, action="store_true", help='Data is nanopore [pacbio]')
  p.add_option('--java', default='/usr/lib/jvm/java-1.8.0/bin/java', help="Path to java 1.8 interpreter [/usr/lib/jvm/java-1.8.0/bin/java]")
  p.add_option('--canu', default='/data/_pipeline_software/canu-1.6/Linux-amd64/bin/canu', help="Path to canu [/data/_pipeline_software/canu-1.6/Linux-amd64/bin/canu]")
  p.add_option('--blast', default='/home/slovett/scripts/blast/split_blasts_multi_faster.py', help="Path to BLAST script [/home/slovett/scripts/blast/split_blasts_multi_faster.py]")
  p.add_option('--trim', default='/data/_scripts/slovett/trim_shift_BAC.py', help="Path to trim shift script [/data/_scripts/slovett/trim_shift_BAC.py]")
  
  opts, args = p.parse_args()

  with open(opts.log, "a") as fout:
    fout.write("###\n%s\n" % " ".join(sys.argv))

  canu_assem = run_canu(opts.input, opts.name, opts.log, opts.java, opts.canu, opts.nanopore)
  
  longest_canu_contig = get_longest_contig(canu_assem)

  backbone_coords = blast_backbone(longest_canu_contig, opts.backbone, opts.log, opts.blast)
  redundancy = find_redundant_ends(longest_canu_contig, opts.log, opts.blast)

  trim_shift(longest_canu_contig, opts.output, backbone_coords, redundancy, opts.log, opts.trim)

  with open(opts.log, "a") as fout:
    fout.write("\n")


def run_canu(reads, assem_name, log, java, canu, nanopore):
  if nanopore:
    platform = "nanopore"
  else:
    platform = "pacbio"
  
  cmd = "%s -p %s -d %s java=%s genomeSize=200k -%s-raw %s maxMemory=1000 maxThreads=100 minMemory=500 minThreads=50" % (canu, assem_name, assem_name, java, platform, reads)
  with open(log, "a") as fout:
    fout.write(cmd+"\n")
  canu=Popen(cmd, shell=True)    
  canu.wait()
  
  return "%s/%s.contigs.fasta" % (assem_name, assem_name)


def get_longest_contig(canu_assem):
  largest_size = 0
  
  with open(canu_assem, "r") as fin:
    for record in SeqIO.parse(fin, "fasta"):
      if len(str(record.seq)) > largest_size:
        largest = (record.id, str(record.seq))
        largest_size = len(str(record.seq))

  with open(canu_assem+".longest", "w") as fout:
    fout.write(">%s\n%s\n" % (largest[0], largest[1]))

  return canu_assem+".longest"
    
  
def blast_backbone(assem_name, backbone_fasta, log, blast):
  cmd = "%s -q %s --ns %s --numHits 10 --numHsps 10" % (blast, backbone_fasta, assem_name)
  with open(log, "a") as fout:
    fout.write(cmd+"\n")
  blast=Popen(cmd, shell=True)    
  blast.wait()
  
  backbone_blast = "%s_blastn_me_%s_parsed.txt" % (backbone_fasta, "/".join(assem_name.split("/")[1:]))
  
  with open(backbone_blast, "r") as fin:
    header = fin.readline()
    found_one_copy = 0
    coords = 0
    for line in fin:
      cols = line.strip().split("\t")
      if float(cols[4])/float(cols[1]) > 0.9:
        if found_one_copy:
          print cols
          print "I already found a copy of the backbone, you need to do this manually"
          sys.exit()
        found_one_copy = 1

        if cols[7] > cols[8]:
          rev = 1
          begin = cols[8]
          end = cols[7]
        else:
          rev = 0
          begin = cols[7]
          end = cols[8]

        coords = (begin, end, rev)

  if not coords:
    print "I did not find the backbone, you need to do this manually"
    sys.exit()

  return coords


def find_redundant_ends(assem_name, log, blast):
  cmd = "%s -q %s --ns %s --numHits 10 --numHsps 10" % (blast, assem_name, assem_name)
  with open(log, "a") as fout:
    fout.write(cmd+"\n")
  blast=Popen(cmd, shell=True)    
  blast.wait()
  
  self_blast = "%s/%s_blastn_me_%s_parsed.txt" % (assem_name.split("/")[0],assem_name.split("/")[1],assem_name.split("/")[1])

  with open(self_blast, "r") as fin:
    header = fin.readline()
    self_hit = fin.readline()

    for line in fin:
      cols = line.strip().split("\t")
      query_len = int(cols[1])
      query_start = int(cols[5])
      query_end = int(cols[6])
      sub_start = int(cols[7])
      sub_end = int(cols[8])
            
      if query_start == 1 and (sub_start == query_len or sub_end == query_len):
        return int(cols[4])

  print "Reached end of self blast without finding redundant ends. You need to do this manually"
  sys.exit()


def trim_shift(assem_name, output, backbone_coords, redundancy, log, trim):
  with open(assem_name, "r") as fin:
    seq_name = fin.readline()[1:]

  breaks_file = assem_name+".breaks"

  with open(breaks_file, "w") as fout:
    fout.write("%s\t%s\t%s\t%s\t%s\n" % (seq_name.strip(), backbone_coords[0], backbone_coords[1], backbone_coords[2], redundancy))

  cmd = "%s -i %s -o %s -b %s" % (trim, assem_name, output, breaks_file)
  with open(log, "a") as fout:
    fout.write(cmd+"\n")
  trim=Popen(cmd, shell=True)    
  trim.wait()


if __name__ == '__main__':
  main()