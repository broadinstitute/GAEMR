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
from gaemr.CoordsFile import CoordsFile
from gaemr.SimpleTable import SimpleTable

parser = OptionParser(usage="usage: %prog [options] nucmer.coords")

parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
                  help='Output table file header name (default=%default)')
parser.add_option('--novel', '-n', action="store_true", default=False, dest="novel",
                  help='Print detailed novel (unaligned) regions (default=%default)')
parser.add_option('--contig_detail', '-c', action="store_true", default=False, dest="details",
                  help='Print detailed contig coverage stats (default=%default)')                  
#parser.add_option('--chart_output', '-t', action="store", default=None, dest="chart_output",
#                  help='Cumulative chart output name (should end with png) (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')


(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply nucmer coords file.")
    sys.exit(-1)

def __get_table_headers(c):
    return ["Fasta File Id", c.get_ref_name(), c.get_query_name()]

def __get_table_title():
    return "Comparison To Reference"

def get_n50(size_list):
    size_list.sort()
    size_list.reverse()
    half = sum(size_list)/2.0
    total = 0
    for i in size_list:
        total += i
        if total >= half:
            return i

    return "NA"
    
def __get_table_data(c):
    table_data = []
    
    table_data.append(["Total Length (bp)", c.get_ref_length(), c.get_query_length()])
    ref_gaps = c.get_ref_gap_sizes()
    query_gaps = c.get_query_gap_sizes()
    table_data.append(["Total Novel Regions", str(len(ref_gaps)), str(len(query_gaps))])
    table_data.append(["Total Novel Bases (bp)", str(sum(ref_gaps)), str(sum(query_gaps))])

    tmp = [["Average Novel Region Size (bp)"],["Largest Novel Region Size (bp)"]]
    if len(ref_gaps):
        tmp[0].append(str("%.0f" % (sum(ref_gaps)/float(len(ref_gaps)))))
        tmp[1].append(str(sorted(ref_gaps)[-1]))
    if len(query_gaps):
        tmp[0].append(str("%.0f" % (sum(query_gaps)/float(len(query_gaps)))))
        tmp[1].append(str(sorted(query_gaps)[-1]))
    table_data += tmp    
    table_data.append(["N50 Novel Region Size (bp)", str(get_n50(ref_gaps)), str(get_n50(query_gaps))])
    table_data.append(["Pct Covered", str(c.get_pct_ref_covered()), str(c.get_pct_query_covered())])
    table_data.append(["Pct Identity", str(c.get_pct_identity()), str(c.get_pct_identity())])
    return table_data

def __get_details_info(base, extension, title, headers, data_function, coord_obj=None):
    output = None
    if base:
        output = base + extension
    if coord_obj:
        data = data_function(coord_obj)
    else:
        data = data_function

    st = SimpleTable(headers,data,title)
    st.print_output(output, options.html)

def main():
    c = CoordsFile(args[0])

    base = None
    if options.output:
        base = options.output

    title = __get_table_title()
    headers = __get_table_headers(c)
    data = __get_table_data(c)

    __get_details_info(base, ".compare_to_ref", title, headers, __get_table_data, c)

    if options.details:
        headers = ["ID", "Length (bp)", "Covered (bp)", "Covered (%)"]
        __get_details_info(base, ".compare_to_ref_query_details", "Query Coverage Details", headers, c.get_query_coverage_details())
        __get_details_info(base, ".compare_to_ref_reference_details", "Reference Coverage Details", headers , c.get_ref_coverage_details())

    if options.novel:
        headers = ["ID", "Start", "End", "Length"]
        __get_details_info(base, ".compare_to_ref_query_novel", "Novel Query Sequence Details", headers, c.get_novel_query_details())
        __get_details_info(base, ".compare_to_ref_reference_novel", "Novel Reference Sequence Details", headers , c.get_novel_ref_details())

    sys.exit(0)

if __name__ == "__main__":
    sys.exit(main())

