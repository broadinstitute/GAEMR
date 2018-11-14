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
from optparse import OptionParser
from gaemr.BamFlagstat import BamFlagstat
from gaemr.BamCoverage import BamCoverage
from gaemr.SimpleTable import SimpleTable

parser = OptionParser(usage="usage: %prog [options] bam_file(s)")

parser.add_option('--output', '-o', action="store", default=None, type='string',dest="output",
    help='Output file header for bam stats (default=%default)')
parser.add_option('--physical_coverage', '-p', action="store_true", default=False, dest="phys_cvg",
    help='Whether to print physical cvg stats as well (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
    help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) < 1:
    parser.error("Must supply at least 1 bam file.")
    sys.exit(-1)

def __get_basename(file):
    tmp = re.split("/",file)
    return re.sub(".bam","",tmp[-1])

def __get_cvg_data(phys_cvg):
    title = "Bam Coverage Stats"
    headers = ["Stat"]
    data = []
    first = 1
    for file in args:
        bf = BamFlagstat(file)
        if bf.is_aligned():
            bc = BamCoverage(file,physical_coverage=phys_cvg)
            base_name = __get_basename(file)
            headers.append(base_name)
            tmp = bc.get_coverage_stats()

            if first:
                data = tmp
                first = 0
                continue

            for i in xrange(len(tmp)):
                data[i].append(tmp[i][1])

    return title,headers,data

def main():

    title,headers,data = __get_cvg_data(options.phys_cvg)

    st = SimpleTable(headers,data,title)
    output = None
    if options.output:
        type = "seq_cvg"
        if options.phys_cvg:
            type = "phys_cvg"
        output = options.output + ".bam_" + type + "_stats"
    st.print_output(output,options.html)

    return 0

if __name__ == "__main__":
    sys.exit(main())
