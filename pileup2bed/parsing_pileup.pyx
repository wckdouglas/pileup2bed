from __future__ import print_function
import string
import sys
from collections import Counter
from functools import partial
import numpy as np
import re
cimport numpy as np
from multiprocessing import Pool
from cpython cimport bool
from itertools import imap

complement_base = string.maketrans('actgnACTGN', 'TGACNTGACN')
np_ord = np.vectorize(ord,otypes=[np.int8])

def strandedBase(str bases_string):
    cdef:
        np.ndarray bases_negative
        np.ndarray bases_positive

    # extract lower and upper strand reads
    bases_negative = np.array(re.findall('[actg]', bases_string))
    bases_positive = np.array(re.findall('[ACTG]', bases_string))
    return bases_positive, bases_negative


cpdef str qualityBases(str bases, str quals, int qual_threshold):
    cdef:
        np.ndarray bases_list
        np.ndarray quals_list
        str out_bases

    # extract high quality bases by numpy array
    bases_list = np.array(list(bases))
    quals_list = np_ord(list(quals)) - 33
    assert len(bases_list) == len(quals_list), 'bases != quals ' + '\n' + \
                                                ''.join(bases) + '!\n' + \
                                                ''.join(quals)
    quality_bases = bases_list[quals_list >= qual_threshold]
    out_bases = ''.join(quality_bases)
    return out_bases


cpdef str printLine(str chrom, int cov_threshold, str pos,
                    bool mismatch_only, str ref,
                    np.ndarray bases, str strand,
                    np.ndarray deletions, np.ndarray insertions):
    cdef:
        str bed_line_str
        int cov

    # print output line as bed format with base count and indel count
    bases_count = Counter(''.join(bases).upper())
    cov = np.sum(bases_count.values())
    if cov >= cov_threshold:
        is_mismatch_only = mismatch_only and cov != bases_count[ref]
        if is_mismatch_only or not mismatch_only:
            bed_line = map(str, [chrom, pos, int(pos)+1,
                                ref, cov, strand,
                                bases_count['A'], bases_count['C'],
                                bases_count['T'], bases_count['G'],
                                len(deletions), len(insertions)])
            bed_line_str = '\t'.join(bed_line)
            return bed_line_str
        else:
            return 'No'
    else:
        return 'No'

def parseBases(str bases, str ref):
    cdef:
        int i = 0
        int j = 0
        int i_max = len(bases)
        int insertion = 0
        int deletion = 0
        str _bases = ''
        str insertion_bases = ''
        str deletion_bases = ''
        str c = ''

    # loop over the bases field
    # if it is ./, add to reference
    # skip for ^ or $ (start and end)
    # insertion/deletion if +/-, follwing number and bases are the indel bases
    # others are mapped base
    while i < i_max:
        c = bases[i]
        if c == '.':
            _bases += ref
        elif c == ',':
            _bases += ref.translate(complement_base).lower()
        elif c == '^':
            i += 1
        elif c == '+':
            j = i + 1
            insert_count = 0
            insert_digit = ''
            while bases[j].isdigit():
                insert_digit += bases[j]
                j += 1
            insert_count = int(insert_digit)
            i = j + insert_count - 1
            insertion_bases += bases[(j):i+1]
            insertion += insert_count
        elif c == '-':
            j = i + 1
            del_count = 0
            del_digit = ''
            while bases[j].isdigit():
                del_digit += bases[j]
                j += 1
            del_count = int(del_digit)
            i = j + del_count - 1
            deletion_bases += bases[(j):i+1]
            deletion += del_count
        elif c != '$':
            _bases += c
        i += 1
    return _bases, insertion, deletion, insertion_bases, deletion_bases


cpdef str processLine(int qual_threshold, int cov_threshold,
                    bool mismatch_only, str line):
    # define variables
    cdef:
        # define split input from mpileup line
        str chrom, pos, ref, cov, inbases, quals
        # define processed line result
        str bases, insertion_bases, deletion_bases
        int coverage, insertion, deletion
        np.ndarray bases_positive, bases_negative
        str print_string

    # split mpileup line
    line = line.strip('\n')
    fields = line.split('\t')
    chrom,  pos, ref, cov, inbases, quals = fields
    coverage = int(cov)

    if coverage > cov_threshold:
        # using bases field to get information
        result = parseBases(inbases, ref)
        bases, insertion, deletion, insertion_bases, deletion_bases = result
        assert len(bases) == len(quals), 'Wrongly parsed!! ' +\
                                         bases + ' ' + \
                                         str(coverage)

        # extract high quality bases only
        bases = qualityBases(bases, quals, qual_threshold)

        # identify upper and lower strand reads
        bases_positive, bases_negative = strandedBase(bases)

        # identify indel from upper and lower strand
        insertion_positive, insertion_negative = strandedBase(insertion_bases)
        deletion_positive, deletion_negative = strandedBase(deletion_bases)

        # print line
        printFunc = partial(printLine, chrom, cov_threshold, pos, mismatch_only)
        ref = ref.upper()
        print_list = map(printFunc, [ref.translate(complement_base), ref],
                                 [bases_negative, bases_positive],
                                 ['-', '+'],
                                 [deletion_negative, deletion_positive],
                                 [insertion_negative, insertion_positive])
        print_string = '\n'.join(print_list).strip('No').strip('\n')
        return print_string
    else:
        return 'No'


cpdef int analyzeFile(str filename, int qual_threshold,
                    int cov_threshold, bool mismatch_only, int threads):
    cdef:
        int lineno
        str line, print_string

    header = ['chrom', 'start', 'end', 'ref_base', 'coverage', 'strand',
              'A', 'C', 'T', 'G', 'deletions', 'insertions']
    print('\t'.join(header), file=sys.stdout)
    lineFunc = partial(processLine, qual_threshold, cov_threshold, mismatch_only)

    if filename == '-':
        handle = sys.stdin
    else:
        handle = open(filename, 'r')

    # runnning process
    pool = Pool(threads)
    processes = pool.imap(lineFunc, handle)
    for lineno, print_string in enumerate(processes):
        if print_string != 'No' and print_string != '':
            print(print_string, file=sys.stdout)
        if lineno % 100000 == 0:
            print('Parsed %i lines' % (lineno), file=sys.stderr)
    pool.close()
    pool.join()

    if filename != '-':
        handle.close()
    return 0
