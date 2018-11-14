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
import re
from gaemr.RunCommand import RunCommand
from optparse import OptionParser
from gaemr.BamFlagstat import BamFlagstat
from gaemr.SimpleTable import SimpleTable

parser = OptionParser(usage="usage: %prog [options] bam_file(s)")

parser.add_option('--output', '-o', action="store", default=None, type='string',dest="output",
                   help='Output file header for bam stats (default=%default)')
parser.add_option('--qc_pass', '-p', action="store_true", default=None, dest="qc_pass",
                   help='Whether to print qc pass reads only (default=%default)')
parser.add_option('--qc_fail', '-f', action="store_true", default=None, dest="qc_fail",
                   help='Whether to print qc fail reads only (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                   help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.error("Must supply at least 1 bam file.")
    sys.exit(-1)

def __check_qc_options(qc_pass,qc_fail):
    if qc_pass:
        if qc_fail:
            return None
        return 1
    if qc_fail:
        return 0
    return None

def __get_header_from_qc(qc_pass):
    if qc_pass is None:
        return "All Reads"
    if qc_pass:
        return "QC-Passed"
    if not qc_pass:
        return "QC-Failed"

def __get_basename(file):
    tmp = re.split("/",file)
    return re.sub(".bam","",tmp[-1])

def __get_table_data(qc_pass):
    title = "Simple Bam Stats"
    headers = ["Stat"]
    data = []
    first = 1
    for file in args:
        b = BamFlagstat(file)
        base_name = __get_basename(file)
        headers.append(base_name + " (" + __get_header_from_qc(qc_pass) + ")")
        if b.is_aligned():
            tmp = b.get_aligned_stats_table(qc_pass)
        else:
            tmp = b.get_unaligned_stats_table(qc_pass)
        if first:
            data = tmp
            first = 0
            continue
        for i in xrange(len(tmp)):
            data[i].append(tmp[i][1])

    return title,headers,data

def main():

    qc_pass = __check_qc_options(options.qc_pass,options.qc_fail)

    title,headers,data = __get_table_data(qc_pass)

    st = SimpleTable(headers,data,title)
    output = None
    if options.output:
        output = options.output + ".simple_bam_stats"
    st.print_output(output,options.html)

    return 0

if __name__ == "__main__":
    sys.exit(main())
