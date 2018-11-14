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
from gaemr.Assembly import Assembly
import gaemr.PlatformConstant as pc

constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option('--output_header', '-o', action="store", default="assembly", type='string',dest="output",
    help='An output header for output files (default=%default)')
parser.add_option('--contigs', '-C', action="store", default=None, type='string',dest="contigs",
    help='Contigs fasta file (default=%default)')
parser.add_option('--scaffolds', '-S', action="store", default=None, type='string',dest="scaffolds",
    help='Scaffolds fasta file (default=%default)')
parser.add_option('--agp', '-A', action="store", default=None, dest="agp",
    help='AGP filename used to build up the assembly (default=%default)')
parser.add_option('--min_contig', '-c', action="store", default=1, type='int', dest="minContig",
    help='Minimum contig length to include in stats (default=%default)')
parser.add_option('--min_scaffold', '-s', action="store", default=1, type='int', dest="minScaffold",
    help='Minimum scaffold length to include in stats (default=%default)')
parser.add_option('--min_gap_size', '-g', action="store", default=constant.MIN_GAP_SIZE, type='int', dest="minGap",
    help='Minimum length of stretch of Ns to consider gap sequence (default=%default)')
parser.add_option('--no_agp_version_2', '-n', action="store_false",default=True, dest="want_agp_version_2",
    help='Do not output agp format in version 2.0 format (default=%default)')
parser.add_option('--rename', '-r', action="store_true",default=False,dest="rename",
    help='Whether or not to rename contig/scaffolds to current NCBI conventions (default=%default)')

(options, args) = parser.parse_args()

if not options.scaffolds and not (options.contigs and options.agp):
    parser.error("Must supply either scaffolds file or a contigs and agp file.")
    sys.exit(1)

def main():
    contigs_output = options.output + ".contigs.fasta"
    scaffold_output = options.output + ".scaffolds.fasta"
    agp_output = options.output + ".agp"

    if options.contigs and options.agp:
        a = Assembly(options.contigs,None,None,options.minGap,options.minContig,options.minScaffold,options.agp)
    else:
        a = Assembly(options.scaffolds,None,None,options.minGap,options.minContig,options.minScaffold,None)

    return_code = a.print_assembly(file=scaffold_output, rename=options.rename)
    if not return_code:
        a.print_assembly_contigs(output=contigs_output, rename=options.rename)
        a.print_assembly_agp(output=agp_output,want_version_2=options.want_agp_version_2,rename=options.rename)

    return 0

if __name__ == "__main__":
    sys.exit(main())
