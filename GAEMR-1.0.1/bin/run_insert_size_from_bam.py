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


parser = OptionParser(usage="usage: %prog [options] bam_file")

parser.add_option('--output_header', '-o', action="store", default=None, type='string',dest="output",
    help='Output header for Picard insert size files (default=%default)')
parser.add_option('--insert_size', '-i', action="store", default=None, type='int', dest="insert_size",
    help='Expected insert size (default=%default)')
parser.add_option('--std_dev', '-s', action="store", default=None, type='int', dest="std_dev",
    help='Expected insert size standard deviation (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply a bam file.")
    sys.exit(-1)

def __get_file_base_name(file):
    if re.search("/",file):
        tmp = file.split("/")
        return tmp[-1]
    return file

def __make_prefix_from_files(file):
    basename = __get_file_base_name(file)
    return re.sub(".bam$","",basename)

def __build_insert_size_command(file,header,insert_size,std_dev):
    insert_size_cmd = constant.PICARD + "CollectInsertSizeMetrics.jar"
    cmd = ["java","-jar",insert_size_cmd,"I=",file]
    if insert_size:
        if std_dev:
            cmd += ["W=",str(insert_size + (3*std_dev))]
        else:
            cmd += ["W=",str(int(insert_size * 2.0))]
    histo = header + ".insert_size.histo.pdf"
    metrics = header + ".insert_size.metrics"
    cmd += ["H=",histo,"O=",metrics,"TMP_DIR=",constant.PICARD_TMP_DIR]
    return cmd

def main():

    output_header = options.output
    if not output_header:
        output_header = __make_prefix_from_files(args[0])

    rc = RunCommand(__build_insert_size_command(args[0],output_header,options.insert_size,options.std_dev))
    out = rc.run_command()

    return 0

if __name__ == "__main__":
    sys.exit(main())
