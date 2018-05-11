#!/usr/bin/env python

from __future__ import division
import optparse
import numpy as np
import gzip

def main():
    #To parse command line
    usage = """usage: %prog [options]
    If processing a single sample, the fastq files can be specified using the -1, -2, and -I flags.
    Alternatively, you can supply an extra file as an argument that should contain one line per sample to process.
    Each line should contain three tab-delimited values: read1.fastq, read2.fastq, index.fastq
    """
    p = optparse.OptionParser(usage)

    p.add_option('-1', '--read1', help='Read 1 fastq file [None]')
    p.add_option('-2', '--read2', help='Read 2 fastq file [None]')
    p.add_option('-I', '--index', help='Index fastq file [None]')
    p.add_option('-q', '--minQual', type='float', default='30', help='Minimum average quality of index for a pair ot be maintained [30]')
    p.add_option('--offset', type='int', default='33', help='ASCII quality score offset [33]')
    opts, args = p.parse_args()
    
    if opts.read1 and opts.read2 and opts.index: filter_index(opts)
    elif args:
        fin=open(args[0], 'r')
        for line in fin:
            opts.read1, opts.read2, opts.index = line.strip().split('\t')
            filter_index(opts)

    else: print 'Must either specify fastq files with -1, -2 and -I flags, or with an extra file that specifies the fastq files for one or more samples'
    
#read1 and read2 should be two different fastq files, containing the two reads from a paired end run
#They should contain only paired seqs and they should be in the same order in the two files
#Each read must only span 4 lines in a fastq
def filter_index(opts):
    
    #Open new output files
    fout_r1 = gzip.open('%s_Q%d.fastq.gz' % (opts.read1.split('.')[0], opts.minQual), 'wb')
    fout_r2 = gzip.open('%s_Q%d.fastq.gz' % (opts.read2.split('.')[0], opts.minQual), 'wb')
    fout_i = gzip.open('%s_Q%d.fastq.gz' % (opts.index.split('.')[0], opts.minQual), 'wb')
    #Open input fastq files
    if opts.read1.endswith("gz"):
      fin_r1 = gzip.open(opts.read1, 'rb')
    else:
      fin_r1 = open(opts.read1, 'r')
    if opts.read2.endswith("gz"):
      fin_r2 = gzip.open(opts.read2, 'rb')
    else:
      fin_r2 = open(opts.read2, 'r')
    if opts.index.endswith("gz"):
      fin_i = gzip.open(opts.index, 'rb')
    else:
      fin_i = open(opts.index, 'r')
    
    good_index_count=0
    while 1:
        first_mate=[]
        second_mate=[]
        index=[]
        #Extract a single sequence from each read file
        for i in range(4):
            first_mate.append(fin_r1.readline())
            second_mate.append(fin_r2.readline())
            index.append(fin_i.readline())
        if first_mate[0]=='': 
            #Close all open files
            fout_r1.close()
            fout_r2.close()
            fout_i.close()
            fin_r1.close()
            fin_r2.close()
            fin_i.close()
            return
        else: 
            if good_qual(index[3], opts.minQual, opts.offset):
                good_index_count+=1
                fout_r1.writelines(first_mate)
                fout_r2.writelines(second_mate)
                fout_i.writelines(index)
    print 'Number of good pairs remaining: %d' % (good_index_count)

def good_qual(index_quals, min_qual, offset):
    index_quals=[ord(x)-offset for x in index_quals.strip()]
    if np.average(index_quals)>=min_qual: return True
    else: return False


###------------->>>

if __name__ == "__main__":
    main()
