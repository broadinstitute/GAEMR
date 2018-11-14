#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import subprocess

import sys
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()
from optparse import OptionParser
from gaemr.BlastFilter import VecScreenFilter

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option('--input', '-i', action="store", default=None, type='string', dest="input",
    help='input file  (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
    help='output file header (default=%default)')
parser.add_option('--threads', '-t', action="store", default="4", type='int', dest="threads",
    help='Number of threads to use for blastn (default=%default)')
parser.add_option('--illumina_primer', '-I', action="store_true", default=False, dest="primer_only",
    help='Only report hits to Illumina adapter/primer (default=%default)')

(options, args) = parser.parse_args()

# check inputs or die
if not options.input:
    parser.error("Must supply input file with --input.")
    sys.exit(-1)

if not options.output:
    parser.error("Must supply output file with --output.")
    sys.exit(-1)


def command(cmd):
    print cmd
    subprocess.check_call(cmd, shell=True)


def main():

    command(constant.GAEMR + "/run_blast.py --vecscreen --output " + options.output + ".vecscreen.xml "+ constant.BLAST_UNIVEC + " " + options.input)
    vec = VecScreenFilter(options.output + ".vecscreen.xml")
    vec.filter_weak_hits()
    if options.primer_only:
        vec.select_Illumina_hits()
    vec.remove_redundant_hits()
    vec_list = vec.get_removal_coords()
    if len(vec_list) == 0:
        print "No Valid hits to vector found."

    final_string = vec.get_m8(0)


    ###Write out contamination
    contam_output = open(options.output + ".vecscreen_filtered.m8", 'w')
    contam_output.write(final_string)
    contam_output.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
