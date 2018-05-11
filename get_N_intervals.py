#!/usr/bin/env python

import optparse, os
from biolad import read_fasta_dict

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-i', '--input', help='Input fasta [None, REQD]')

    opts, args = p.parse_args()

    get_N_intervals(opts)


def get_N_intervals(opts):
    input_dict = read_fasta_dict(opts.input)

    for name, seq in input_dict.iteritems():
        print "%s\nStart\tEnd\tSize" % (name)
        last_was_n = 0
        for i in xrange(len(seq)):
            if seq[i] is 'N':
                if last_was_n == 0:
                    range_start = i+1
                last_was_n = 1
            elif last_was_n == 1:
                print "%s\t%s\t%s" % (range_start, i, i-range_start+1)
                last_was_n = 0


###---------->>>

if __name__ == '__main__':
    main()
