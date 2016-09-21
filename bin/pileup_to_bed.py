#!/usr/bin/env python

from __future__ import print_function
import sys
import string
import argparse
import pyximport
import numpy as np
from pileup2bed.parsing_pileup import analyzeFile
pyximport.install(setup_args={'include_dirs': np.get_include()})


def getopt():
    descriptions = 'This program tries to convert samtools mipleup' + \
                ' format to base count table with' + \
                ' coverage threshold and quality threshold'
    parser = argparse.ArgumentParser(description=descriptions)
    parser.add_argument('-i', '--input', required=True,
                        help='mpileup format input filename (default: <->)')
    parser.add_argument('-c', '--cov', default=10, type=int,
                        help='coverage threshold (default: 10)')
    parser.add_argument('-q', '--qual', default=33, type=int,
                        help='quality threshold (default: 33)')
    parser.add_argument('-m','--mismatch_only', action='store_true', help='output position with mismatch only')
    parser.add_argument('-t','--threads', default=1, help='threads to use (default = 1)', type=int)
    args = parser.parse_args()
    return args


def main():
    # get arguments and call cython function
    args = getopt()
    qual_threshold = args.qual
    cov_threshold = args.cov
    filename = args.input
    mismatch_only = args.mismatch_only
    threads = args.threads
    header = ['chrom', 'start', 'end', 'ref_base', 'coverage', 'strand',
              'A', 'C', 'T', 'G', 'deletions', 'insertions']
    print('\t'.join(header), file=sys.stdout)
    handle = sys.stdin if filename == '-' else open(filename, 'r')

    analyzeFile(handle, qual_threshold, cov_threshold, mismatch_only, threads)

if __name__ == '__main__':
    main()
