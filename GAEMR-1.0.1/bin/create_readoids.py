#! /usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser


parser = OptionParser(usage="usage: %prog [options] <fasta file>")

parser.add_option('--length', '-l', action="store", default="100", type='int',dest="rd_ln",
                  help='Readoid length (default=%default)')
parser.add_option('--distance', '-d', action="store", default="100000", type='int',dest="distance",
                  help='Distance between readoids (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
                  help='Output fasta file (default=%default)')
parser.add_option('--window_size', '-w', action="store", default="10000", type='int', dest="window_size",
                  help='Window size to slide across fasta file (default=%default)')

(options, args) = parser.parse_args()

# check inputs or die
if len(args) !=1:
    parser.error("Must supply fasta file.")
    sys.exit(-1)

def makeReadoids(input,output):
    try: 				# open up output
        out_fh = open(output, 'w')
        in_fh  = open(input, 'r')
        pair_count = 0
        for record in SeqIO.parse(in_fh,'fasta'):				# we are telling module it is fasta
            seq = record.seq
            id  = record.id
            ln  = len(seq)
            last_position = ln - options.distance
            rd_ln = options.rd_ln
            distance = options.distance
            window_size = options.window_size

            for i in range(0,last_position,window_size):
                read1     = seq[i:i + rd_ln]
                rd2_start = i + (distance - rd_ln)
                rd2_end   = rd2_start + rd_ln
                read2     = seq[rd2_start:rd2_start + rd_ln]
                read2_rc  = read2.reverse_complement()
                header1   = 'p' + str(pair_count) + '_' + id + '_' + str(i)	 # cant combine strings and ints, hence use of str()
                header2   = 'p' + str(pair_count) + '_' + id + '_' + str(rd2_end) # cant combine strings and ints, hence use of str()
                out_fh.write(">%s\n%s\n" % (header1, read1))
                out_fh.write(">%s\n%s\n" % (header2, read2_rc))
                pair_count += 1




    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1

    # see ~/dev/ScaffoldAssessment.pl makeReadoids perl sub

def main():

    output_fasta=options.output
    if not output_fasta:
        output_fasta=re.sub(".fasta",".readoid.fasta",args[0])

    makeReadoids(args[0],output_fasta)

if __name__ == "__main__":			# this will call main and start the script
    sys.exit(main())

