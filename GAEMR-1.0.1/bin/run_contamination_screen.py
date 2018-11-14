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
from gaemr.BlastFilter import VecScreenFilter, gContamFilter, MitoFilter, rRNAFilter, CombinedContaminationObject, BlastFilter

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option('--input', '-i', action="store", default=None, type='string', dest="input",
    help='Input file fasta')

parser.add_option('--output', '-o', action="store", default="output", type='string', dest="output",
    help='Output file header: ')

parser.add_option('--output_format', '-f', action="store", default="m8", type='string', dest="output_format",
    help='Output file format: (default=%default)\n\tm8 - BLAST m8 format\n\tSubPrep - SubmissionPrep format contig_id, contig length, contam start, contam end, hit description\n\tgeneral - General removal format contig_name, contig length, contam start, contam end')

parser.add_option('--stdout', '', action="store_true", default=False, dest="stdout",
    help='Output to stdout (default= False)')

parser.add_option('--redundant', '-r', action="store_true", default=False, dest="redundant",
    help='Remove redundant hits (default=%default)')

parser.add_option('--extend', '-e', action="store_true", default=False, dest="extend",
    help='Merge and extend overlapping hits (default=%default)')

parser.add_option('--reach', '-R', action="store", type='int', default=200, dest="reach",
    help='Size of reach (distance between hits) for extension and merging (default=%default)')

parser.add_option('--threads', '-t', action="store", default="4", type='int', dest="threads",
    help='Number of threads to use for blastn (default=%default)')

### All
parser.add_option('--all', '', action="store_true", default=False, dest="all",
    help='Run all screens: gcontam, mitochondria, rRNA and Vecscreen (default=%default)')


### gContam
parser.add_option('--gcontam', '-g', action="store_true", default=False, dest="gcontam",
    help='Run Screen against NCBI\'s gcontam db (default=%default)')


### Mito
parser.add_option('--mito', '-m', action="store_true", default=False, dest="mito",
    help='Run mitochondria screen (default=%default)')

parser.add_option('--pct_mito', '-p', action="store", default=".30", type='float', dest="pct_mito",
    help='Fraction of contig required to hit mitochondria to be filtered [default: .30')

### VecScreen
parser.add_option('--vecscreen', '-v', action="store_true", default=False, dest="vecscreen",
    help='Run NCBI VecScreen (default=%default)')

parser.add_option('--illumina_primer', '-I', action="store_true", default=False, dest="primer_only",
    help='Only report hits to Illumina adapter/primer in VecScreen (default=%default)')

### rRNA
parser.add_option('--rRNA', '', action="store_true", default=False, dest="rRNA",
    help='Run NCBI rRNA Screen (default=%default)')

###Custom
parser.add_option('--custom', '-c', action="store", default=None, type="string", dest="custom_db",
    help='Custom db for contamination screen')

parser.add_option('--filter_as', '-F', action="store", default=None, type="string", dest="filter_as",
    help='Filter hits from custom database like:\tgcontam\tmitochondria\tvecscreen\trRNA\tNone (default=%default)')

(options, args) = parser.parse_args()

# check inputs or die
if not options.input:
    parser.error("Must supply input file with --input.")
    sys.exit(-1)

if not options.output:
    parser.error("Must supply output file with --output.")
    sys.exit(-1)

if not options.mito and not options.gcontam and not options.vecscreen and not options.all and not options.custom_db and not options.rRNA:
    parser.error("Must select a screen: --all, -g, -m, -v, --rRNA or -c")
    sys.exit(-1)

if options.all:
    options.mito = True
    options.gcontam = True
    options.vecscreen = True
    options.rRNA = True

### Change pct_mito from percent to fraction
if options.pct_mito > 1:
    pct_mito = options.pct_mito / 100
else:
    pct_mito = options.pct_mito

def log_filter(filter, Name):
    hit_list = filter.get_removal_coords()
    print "Filter: " + Name + "\t" + str(len(hit_list))


def command(cmd):
    print cmd
    subprocess.check_call(cmd, shell=True)


def vecscreen_filters(xml):
    vec = VecScreenFilter(xml)
    vec.filter_weak_hits()
    if options.primer_only:
        vec.select_Illumina_hits()
    log_filter(vec, "Vector Screen")
    return vec


def rRNA_filters(xml):
    rRNA = rRNAFilter(xml)
    rRNA.filter_by_length(100)
    rRNA.filter_by_pct_covered(0.8)
    log_filter(rRNA, "rRNA Screen")
    return rRNA


def gcontam_filters(xml):
    gcontam = gContamFilter(xml)
    gcontam.ignore_mito_hits()
    gcontam.filter_hits()
    log_filter(gcontam, "gContam Screen")
    return gcontam


def mito_filters(xml):
    mito = MitoFilter(xml, pct_mito)
    mito.select_valid_mito_hits()
    log_filter(mito, "Mitochondria Screen")
    return mito


def no_filter(xml):
    raw = BlastFilter(xml)
    log_filter(raw, "Raw Hits - No Filter")
    return raw


def list_to_string(list):
    out_string = ""
    for hit in list:
        out_string += (''.join([str(i) + "\t" for i in hit]))
        out_string = out_string[:-1] + "\n"

    return out_string


def main():
    obj_list = []
    if options.mito or options.gcontam:
        command("run_blast.py --ncbi -t " + str(options.threads) + " --output " + options.output + ".ncbi_screen.xml " +
                constant.BLAST_MITOGCONTAM + " " + options.input)
        ## Screen for mito
        if options.mito:
            mito = mito_filters(options.output + ".ncbi_screen.xml")
            obj_list.append(mito)
            ## Screen for gContam
        if options.gcontam:
            gcontam = gcontam_filters(options.output + ".ncbi_screen.xml")
            obj_list.append(gcontam)

    ## Screen for Vector
    if options.vecscreen:
        command("run_blast.py --vecscreen -t " + str(options.threads) + " --output " + options.output + ".vecscreen.xml " +
                constant.BLAST_UNIVEC + " " + options.input)
        vec = vecscreen_filters(options.output + ".vecscreen.xml")
        obj_list.append(vec)

    ## Screen for rRNA
    if options.rRNA:
        command("run_blast.py --rRNA -t " + str(options.threads) + " --output " + options.output + ".rRNAscreen.xml " +
                constant.BLAST_rRNA + " " + options.input)
        rRNA = rRNA_filters(options.output + ".rRNAscreen.xml")
        obj_list.append(rRNA)

    ## Screen custom db
    if options.custom_db:
        command("run_blast.py -b megablast -t " + str(options.threads) + " -d --output " + options.output +
                ".custom_screen.xml " + options.custom_db + " " + options.input)
        if options.filter_as == "gcontam":
            custom = gcontam_filters(options.output + ".custom_screen.xml")
        elif options.filter_as == "mito":
            custom = mito_filters(options.output + ".custom_screen.xml")
        elif options.filter_as == "vecscreen":
            custom = vecscreen_filters(options.output + ".custom_screen.xml")
        elif options.filter_as == "rRNA":
            custom = rRNA_filters(options.output + ".custom_screen.xml")
        elif not options.filter_as:
            custom = no_filter(options.output + ".custom_screen.xml")
        else:
            print "Error filter_as option: " + options.filter_as + " is not an accepted value."
            sys.exit(-1)
        obj_list.append(custom)


    ## Combine hits, remove redundant hits, merge and extend hits
    if len(obj_list) > 1:
        combined = CombinedContaminationObject(obj_list)
    else:
        combined = obj_list[0]

    if options.redundant:
        combined.remove_redundant_hits()
        log_filter(combined, "Post Redundancy")

    if options.extend:
        combined.merge_and_extend_hits(options.reach)
        log_filter(combined, "Post Merging and Extension")

    if options.output_format.lower() == "m8":
        final_string = combined.get_m8()
    elif options.output_format.lower() == "subprep":
        final_string = list_to_string(combined.get_removal_coords(True))
    elif options.output_format.lower() == "general":
        final_string = list_to_string(combined.get_removal_coords(False))
    else:
        print "Error: Output Format can not be None"
        sys.exit(-1)

    ###Write out contamination
    if options.stdout:
        print "HIT LIST:\n" + final_string
    else:
        contam_output = open(options.output + ".contam.list", 'w')
        contam_output.write(final_string)
        contam_output.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
