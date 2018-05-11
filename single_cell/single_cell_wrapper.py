#!/usr/bin/python

'''
Sean Lovett
9/11/2017

This script is meant to simpify running the Wafergen single cell analysis workflow. It will generate config and, if necessary, metadata files and kick off the pipeline.

Minimum input requirements:
Read1 and read2 fastqs (may be gzipped)
Well list OR metadata file
Indexed reference genome
Reference genome GTF

Note that, for the time being, doing analysis on anything except human or mouse requires extra steps.
'''

from __future__ import division
from subprocess import Popen, PIPE
import optparse
import os.path

def main():
  usage = '%prog [options]\nIn case of strange output, please sanity check config and metadata files.'
  p = optparse.OptionParser()

  p.add_option('-n', '--name', help='Sample name [None,REQD]')
  p.add_option('-d', '--desc', default = 'No description given', help='Sample description [NULL]')
  p.add_option('-m', '--meta', help='Metadata file, will be created if a well list is provided [None,REQD]')
  p.add_option('-w', '--wells', help='Well list, required if metadata file not provided [None]')
  p.add_option('-1', '--read1', help='Input R1 [None, REQD]')
  p.add_option('-2', '--read2', help='Input R2 [None, REQD]')
  p.add_option('-o', '--output', default=".",help='Analysis output directory [.]')
  p.add_option('-t', '--threads', type='int', default= 10, help='Number of threads, should generally be left alone [10]')
  p.add_option('-r', '--ref', help='Reference genome index prefix [None, REQD]')
  p.add_option('-g', '--gtf', help='Reference genome gtf [None,REQD]')

  opts, args = p.parse_args()

  # Absolute file paths are required, so we might as well enforce that now.
  opts.read1 = os.path.abspath(opts.read1)
  opts.read2 = os.path.abspath(opts.read2)
  opts.meta = os.path.abspath(opts.meta)
  opts.wells = os.path.abspath(opts.wells)
  opts.output = os.path.abspath(opts.output)
  opts.ref = os.path.abspath(opts.ref)
  opts.gtf = os.path.abspath(opts.gtf)

  # Light sanity checking of user input before we do anything time-consuming
  # At least one of the well list and metadata file must exist, and we don't want to overwrite existing metadata
  if not opts.wells and not os.path.isfile(opts.meta):
    print("Well list is required if metadata file does not exist")
    exit()
  elif opts.wells and os.path.isfile(opts.meta):
    print("Will not overwrite existing metadata file: %s" % opts.meta)
    exit()

  # The pipeline requires unzipped fastqs for some reason
  if opts.read1.endswith("gz"):
    opts.read1 = run_unpigz(opts.read1)
  if opts.read2.endswith("gz"):
    opts.read2 = run_unpigz(opts.read2)

  # Create metadata and config files as necessary
  if opts.wells and not os.path.isfile(opts.meta):
    write_meta(opts)
  config = write_config(opts)
  
  # Run the pipeline and hope you did everything right.
  run_sca(config)


# Run the pipeline. This method is completely self-explanatory.
def run_sca(config):
  cmd = "/data/_pipeline_software/wafergen_pipeline/versions/v1.0.2/bin/sc_wta_analysis.pl -c %s" % config
  sca_run=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
  sca_run.wait()


# Convert the well list to a metadata file. This may need to change at some point, but this converter would have worked for everything so far.
def write_meta(opts):
  meta_lines = []
  with open(opts.wells, "r") as fin:
    header = fin.readline()
    for line in fin:
      cols = "_".join(line.split(" ")).strip().split("\t")
      meta_line = "%s\tN/A\t%s\tMini_Single\tN/A\tN/A\t%s\t%s\tN/A\n" % (cols[5], cols[4], cols[0], cols[1] )
      meta_lines.append(meta_line)
  with open(opts.meta, "w") as fout:
    for line in meta_lines:
      fout.write(line)


def run_unpigz(gzip_file):
  cmd = "unpigz %s" % gzip_file
  new_name = ".".join(gzip_file.split('.')[:-1])
  
  unpigz_run=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
  unpigz_run.wait()

  return new_name


# Writes a minimal config file based on user options. This will need to be sanity checked from time to time.
def write_config(opts):
  out_file = opts.name + ".config"
  
  with open(out_file,"w") as fout:
    fout.write("sample_name\t%s\n" % opts.name)
    fout.write("sample_description\t%s\n" % opts.desc)
    fout.write("barcode_length_1\t11\n")
    fout.write("barcode_length_2\t0\n")
    fout.write("umi_length_1\t10\n")
    fout.write("umi_length_2\t0\n")
    fout.write("metadata_file\t%s\n" % opts.meta)
    fout.write("src_fastq1_file\t%s\n" % opts.read1)
    fout.write("src_fastq2_file\t%s\n" % opts.read2)
    fout.write("aligner\tstar\n")
    fout.write("analysis_dir\t%s\n" % opts.output)
    fout.write("max_concurrent_jobs\t%s\n" % opts.threads)
    fout.write("ref_genome_prefix\t%s\n" % opts.ref)
    fout.write("ref_gtf_file\t%s\n" % opts.gtf)
  return out_file


if __name__=="__main__":
  main()