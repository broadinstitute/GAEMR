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
from gaemr.RunCommand import RunCommand
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options] scaffolds.fasta OR contigs.fasta")

parser.add_option('--classify', '-c', action="store_true", default=None,dest="classify",
                  help='Whether or not to run classifier to find out taxonomy of rna hits (default=%default)')
parser.add_option('--rnammer_output', '-f', action="store", default=None, type='string',dest="rnammer_out",
                  help='Rnammer output filename -- will be fasta format (default=%default)')
parser.add_option('--rdp_classifier_output', '-r', action="store", default=None, type='string',dest="rdp_out",
                  help='RDP Classifier output filename (default=%default)')
parser.add_option('--gene', '-g', action="store", default="16s", type='string', dest="gene",
                  help='Gene of interest: 5s, 16s, 23s or any combo.  Choose "all" for all three. (default=%default)')
parser.add_option('--superkingdom', '-s', action="store", default="bac", type='string', dest="superkingdom",
                  help='Superkingdom of interest: arc, bac, euk or any combo.  Choose "all" for all three. (default=%default)')


(options, args) = parser.parse_args()

if len(args) != 1 or not options.rnammer_out:
    parser.error("Must supply input fasta file and specify rnammer output file.")
    sys.exit(1)

# private function to make sure user put in correct args
def __check_inputs(classify, rdp_out):
    if classify and not rdp_out or rdp_out and not classify:
        return 0
    return 1

# private function to print out command line input error
def __print_parse_error(type, input, argument):
    print "Unrecognized " + type + " name, " + input + ".\nPlease use " + argument + " or any combo separated by comma.\n"
    sys.exit(-1)    

# private function to build gene arg in rnammer
def __get_gene_string(gene):
    gene_string = ""
    data = gene.split(',')

    for i in data:
        if i == '5s':
            gene_string += "tsu" + ","
        if i == '16s':
            gene_string += "ssu" + ","
        if i == '23s':
            gene_string += "lsu" + ","
        if i == 'all':
            gene_string += "tsu,ssu,lsu,"

    if gene_string:
        return gene_string[:-1]
    else:
        __print_parse_error("gene", gene, "5s, 16s, 23s")

# private function to make rnammer superkingdom arg
def __get_superkingdom(superkingdom):
    super_string = ""
    data = superkingdom.split(',')

    for i in data:
        if i == 'all':
            super_string += "arc,bac,euk,"
        if i == 'arc' or i == 'bac' or i == 'euk':
            super_string += i + ","

    if super_string:
        return super_string[:-1]
    else:
        __print_parse_error("superkingdom", superkingdom, "arc, bac, euk")

# private functions to build commands
def __build_rnammer_cmd(input, output, gene, superkingdom):
    return [constant.RNAMMER, "-S", __get_superkingdom(superkingdom), "-m", __get_gene_string(gene), "-f", output, "-multi", input]

def __build_rdp_cmd(input, output):
    return ["java", "-Xmx400m", "-jar", constant.RDP, "-q", input, "-o", output]


def main():

    if __check_inputs(options.classify, options.rdp_out):
        rnammer_cmd = __build_rnammer_cmd(args[0], options.rnammer_out, options.gene, options.superkingdom)
        rc = RunCommand(rnammer_cmd)
        print "RUNNING:  " + rc.get_command()
        rc.run_command()

        if options.classify:
            rdp_cmd = __build_rdp_cmd(options.rnammer_out,options.rdp_out)
            rc = RunCommand(rdp_cmd)
            print "RUNNING:  " + rc.get_command()
            rc.run_command()

        return 0

    else:
        print "If classifying hits, must supply classify flag and rdp output file."
        sys.exit(-1)

if __name__ == "__main__":
    sys.exit(main())


