#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import sys
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()
from gaemr.BamFlagstat import BamFlagstat
from gaemr.Alignment import *

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] -i unmapped.bam -r reference.fasta")

parser.add_option('--input_file', '-i', action="store", default=None, type='string', dest="input",
                  help='Input unaligned BAM/SAM file')

parser.add_option('--output_file_header', '-o', action="store", default=None, type='string', dest="output_header",
                  help='Output file name header [default: input_header.aligned]')

parser.add_option('--reference', '-r', action="store", default=None, type='string', dest="reference",
                  help='reference fasta file')

parser.add_option('--aligner', '-a', action="store", default="BWA", type='string', dest="aligner",
                  help='Aligner {BWA, BowTie} [default: BWA]')

parser.add_option('--short', '-s', action="store_true", default=True, dest="short",
                  help='Align as short reads (<200) [default]')

parser.add_option('--long', '-l', action="store_false", default=True, dest="short",
                  help='Align as long reads (>200)')

parser.add_option('--no_index', '-x', action="store_false", default=True, dest="index",
                  help='Skip building the reference index')

parser.add_option('--mark_duplicates', '-d', action="store_true", default=False, dest="mark_dups",
                  help='Mark duplicates in final aligned bam')

parser.add_option('--threads', '-t', action="store", default=2, type='int', dest="threads",
                  help='Number of threads to be used [default: 2]')

parser.add_option('--tmp_dir', '-T', action="store", default=".", type='string', dest="tmp_dir",
                  help='Tmp dir for alignment files [default: cwd]')

parser.add_option('--unpaired', '-u', action="store_true", default=False, dest="unpaired",
                  help='Force alignment as unpaired [default: False]')


(options, args) = parser.parse_args()

if not options.input:
    parser.error("Must supply an unaligned BAM as input (-i)")
    sys.exit(-1)

if not options.reference:
    parser.error("Must supply a reference (-r)")
    sys.exit(-1)

if not options.output_header:
    options.output_header = options.input.rstrip(".bam").rstrip(".unmapped") + ".aligned"

def main():
    print "INPUT: " + options.input
    print "OUTPUT: " + options.output_header
    print "Aligner: " + options.aligner
    print "Short Reads: " + str(options.short)

    ### Check is unmapped
    #aligned = 0
    bam_stats = BamFlagstat(options.input)
    aligned = bam_stats.is_aligned()
    print "Aligned: " + str(aligned)
    if aligned:
        parser.error("Must supply an UNALIGNED BAM as input (-i) - Please revert and rerun")
        sys.exit(-1)

    #### Check if paired
    #paired = 1
    forced = ""
    if not options.unpaired:
        paired = bam_stats.is_paired()
    else:
        forced = "Forced"
        paired = 0

    print "PAIRED: " + str(paired) + " " + forced

    aligner = options.aligner.upper()
    if aligner == "BWA":
        bwa = BwaAlignment(options.input, options.reference, options.output_header, paired, options.short, threads=options.threads, forced=options.unpaired )
        bwa.run_alignment(options.index, options.tmp_dir, options.mark_dups)
        bwa.clean_up_files()
    elif aligner =="BOWTIE":
        bowtie = BowTieAlignment(options.input, options.reference, options.output_header, paired, options.short, options.threads)
        bowtie.run_alignment(options.index,options.tmp_dir)
        bowtie.clean_up_files()

    return 0

if __name__ == "__main__":
    sys.exit(main())

