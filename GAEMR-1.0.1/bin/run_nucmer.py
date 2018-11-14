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
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()


parser = OptionParser(usage="usage: %prog [options] reference_file query_file")

parser.add_option('--prefix', '-p', action="store", default=None, type='string',dest="prefix",
                   help='Prefix for nucmer output (default=%default)')
parser.add_option('--uniqueness', '-u', action="store", default="mum", type='string', dest="unique",
                  help='Type of anchoring:  mum, mumcand, mumreference (default=%default)')
parser.add_option('--maxmatch', '-x', action="store_true", default=False, dest="maxmatch",
                  help='Run nucmer using -maxmatch option (default=%default)')
parser.add_option('--minmatch', '-i', action="store", type=int, default=None, dest="minmatch",
                  help='Run nucmer using -minmatch <int> option (default=%default)')
parser.add_option('--nosimplify', '-n', action="store_true", default=False, dest="nosimplify",
                  help='Run nucmer using --nosimplify to look for repeats (default=%default)')
parser.add_option('--mummerplot', '-m', action="store_true", default=False, dest="mummerplot",
                  help='Run mummerplot on resulting nucmer output (default=%default)')
parser.add_option('--nofilter', '-f', action="store_false", default=True, dest="filter",
                  help='Only display best alignments in mummerplot (default=%default)')
parser.add_option('--nolayout', '-l', action="store_false", default=True, dest="layout",
                  help='Print mummerplot in logical layout (default=%default)')
parser.add_option('--postscript', '-s', action="store_true", default=False, dest="postscript",
                  help='Whether or not to generate a postscript file in mummerplot (default=%default)')
parser.add_option('--promer', action="store_true", default=False, dest="promer",
                  help='Run promer to do alignment rather than nucmer (default=%default)')
(options, args) = parser.parse_args()

if len(args) != 2:
    parser.error("Must supply reference and query fasta files.")
    sys.exit(-1)

def __get_file_base_name(file):
    if re.search("/",file):
        tmp = file.split("/")
        return tmp[-1]
    return file

def __make_prefix_from_files(ref,query):
    ref_pre = re.sub(".fasta","",__get_file_base_name(ref))
    query_pre = re.sub(".fasta","",__get_file_base_name(query))
    return ref_pre + "_vs_" + query_pre

def __correct_uniq_arg(arg):
    return arg == 'mum' or arg == 'mumcand' or arg == 'mumreference'

def __get_args(uniq,max,min,simplify):
    arg_list = []
    if max:
        arg_list += ["--maxmatch"]
    else:
        if __correct_uniq_arg(uniq):
            arg_list += ["-" + uniq]
    if min:
        arg_list += ["--minmatch",int(min)]
    if simplify:
        arg_list += ["--nosimplify"]
    return arg_list

def __build_nucmer_command(ref,query,prefix,arg_list):
    program = constant.NUCMER
    if options.promer:
        program = constant.PROMER
    return [program,"-p",prefix] + arg_list + ["-o",ref,query]

def __build_mummerplot_command(delta,arg_list,out_type):
    cmd_list = [constant.MUMMERPLOT] + arg_list + [delta]
    if out_type:
        cmd_list += [out_type]
    return cmd_list

def __get_plot_args(filter,layout,prefix):
    arg_list = []
    if filter:
        arg_list += ["-f"]
    if layout:
        arg_list += ["-l"]
    if prefix:
        arg_list += ["-p",prefix]
    return arg_list

def main():

    prefix = options.prefix
    if not options.prefix:
        prefix = __make_prefix_from_files(args[0],args[1])
    arg_list = __get_args(options.unique,options.maxmatch,options.minmatch,options.nosimplify)
    rc = RunCommand(__build_nucmer_command(args[0],args[1],prefix,arg_list))
    out = rc.run_command()

    if options.mummerplot:
        delta_file = prefix + ".delta"
        plot_arg_list = __get_plot_args(options.filter,options.layout,prefix)
        out_type = '--png'
        if options.postscript:
            out_type = '--postscript'
        rc = RunCommand(__build_mummerplot_command(delta_file,plot_arg_list,out_type))
        out = rc.run_command()

    return 0

if __name__ == "__main__":
    sys.exit(main())
