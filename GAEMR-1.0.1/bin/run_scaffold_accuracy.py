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
from gaemr.RunCommand import RunCommand
from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] <ref fasta file> <scaffolds fasta file>")

parser.add_option('--length', '-l', action="store", default="100", type='int',dest="rd_ln",
    help='Readoid length (default=%default)')
parser.add_option('--distance', '-d', action="store", default="100000", type='int',dest="distance",
    help='Distance between readoids (default=%default)')
parser.add_option('--window_size', '-w', action="store", default="10000", type='int', dest="window_size",
    help='Window size to slide across fasta file (default=%default)')
parser.add_option('--error', '-e', action="store", default=".05", type='float',dest="error",
    help='Allowable deviation of distance between readoids (default=%default)')
parser.add_option('--output', '-o', action="store", default="ScaffoldAssessment", type='string', dest="output",
    help='Output prefix (default=%default)')
parser.add_option('--table', '-t', action="store", default="ScaffoldAssessment", type='string', dest="table",
    help='Table output file prefix (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 2:
    parser.error("Must supply reference fasta file and scaffolds fasta file.")
    sys.exit(1)

def main():

    readoids = options.output + ".readoids.fasta"
    rc = RunCommand(['create_readoids.py','-l',str(options.rd_ln),'-d',str(options.distance),'-w',str(options.window_size),
                     '-o',readoids,args[1]])
    sys.stderr.write("Running:  " + rc.get_command())
    rc.run_command()

    detail_output = options.output + ".scaffold_accuracy"
    rc = RunCommand(['run_nucmer.py','-p',detail_output, args[0], readoids])
    sys.stderr.write("Running:  " + rc.get_command())
    rc.run_command()

    coords_file = detail_output + ".coords"
    table_output = options.table + ".scaffold_accuracy"
    rc = RunCommand(['assess_readoids_coords.py','-c',coords_file,'-p',readoids,'-f',args[1],'-r',args[0],'-v', detail_output,
                     '-l',str(options.rd_ln),'-d',str(options.distance),'-e',str(options.error),'-o',table_output])
    sys.stderr.write("Running:  " + rc.get_command())
    rc.run_command()

    return 0

if __name__ == "__main__":
    sys.exit(main())
