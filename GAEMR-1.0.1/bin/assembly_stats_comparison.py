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
from optparse import OptionParser
from gaemr.Assembly import Assembly
from gaemr.AssemblyStatsUtil import AssemblyStatsUtil
from gaemr.SimpleTable import SimpleTable
import gaemr.PlatformConstant as pc

constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options] list")

parser.add_option('--min_contig', '-c', action="store", default=constant.MIN_CONTIG_SIZE, type='int', dest="minContig",
                  help='Minimum contig length to include in stats (default=%default)')
parser.add_option('--min_scaffold', '-s', action="store", default=constant.MIN_SCAFFOLD_SIZE, type='int', dest="minScaffold",
                  help='Minimum scaffold length to include in stats (default=%default)')
parser.add_option('--min_gap_size', '-g', action="store", default=constant.MIN_GAP_SIZE, type='int', dest="minGap",
                  help='Minimum length of stretch of Ns to consider gap sequence (default=%default)')
parser.add_option('--output_file', '-o', action="store", default=None, type='string', dest="output",
                  help='Output file name (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')\

(options, args) = parser.parse_args()


# private function to load assembly data for cross-comparison
def __load_assemblies(data):
    assemblies = []
    for i in data:
        assert len(i) == 3,"Invalid line length:  must have 3 comma-seperated records"
        assemblies.append(Assembly(fasta=i[0], name=i[1], assembler=i[2], minGapSize=options.minGap, minConSize=options.minContig, minScaffSize=options.minScaffold))
    return assemblies

def main():
    data = []

    if len(args) != 1:
        parser.error("Must supply assembly list file.")
        sys.exit(1)    
    try: # open file or fail
        f=open(args[0])
        for lines in f.readlines():
            data.append(lines.rstrip('\n').split(','))
        f.close()

        # load assemblies and load stats module
        assemblies = __load_assemblies(data)
        s = AssemblyStatsUtil(assemblies)

        #print s.get_stats()


        title = "Basic Assembly Stats"
        headers = s.get_assembly_names()
        data = s.get_stats()
        st = SimpleTable(headers,data,title)
        output = None
        if options.output:
            output = options.output + ".assembly_stats_comparison"
        st.print_output(output,options.html)

    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1

    return 0

if __name__ == "__main__":
    sys.exit(main())
