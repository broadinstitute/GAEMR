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
from gaemr.BamFile import BamFile
from gaemr.SffFile import SffFile
from gaemr.FastqFile import FastqFile
from gaemr.FastaFile import FastaFile

parser = OptionParser(usage="usage: %prog [options] reads.file")

parser.add_option('--output_bam', '-o', action="store", default=None, type='string',dest="output",
                  help='An output bam filename (default=%default)')
parser.add_option('--default_qual', '-q', action="store", default='20', type='string',dest="quality_score",
                  help='A default quality to give fasta files with no associated quals (default=%default)')
parser.add_option('--direction', '-d', action="store", default=None, type='string',dest="direction",
                  help='The default direction of the input reads if paired (default=%default)')
parser.add_option('--sample', '-s', action="store", default='sample', type='string', dest="sample",
                  help='The sample name.  Used for picard FastqToSam.jar software (default=%default)')
parser.add_option('--quality_format', '-f', action="store", default='Standard', dest="qual_format",
                  help='Quality format for fastq.  Standard fastq, Illumina fastq (+33) or Solexa fastq (+64) (default=%default)')
parser.add_option('--java_heap', action="store", type='string', default='32G', dest="java_heap",
                  help='Java heap size for Picard modules:  examples 2G, 8G, etc (default=%default)')
(options, args) = parser.parse_args()

# check args
if len(args) < 1:
    parser.error("Must supply one or more reads file.")
    sys.exit(1)

# set up our initialize files dict
def __initialize_files():
    files = {}
    files['FASTA'] = []
    files['QUAL'] = []
    files['FASTQ'] = []
    files['SAM'] = []
    files['BAM'] = []
    files['SFF'] = []
    return files

# private function to assign file to a category
def __read_args(args,files):
    for i in args:
        extension = re.sub(".*\.","",i).upper()
        tmp = extension[0] + extension[-1]
        if extension != 'FASTA' and tmp == 'FA':
            extension = 'FASTA'
        if extension != 'FASTQ' and tmp == 'FQ':
            extension = 'FASTQ'
        if extension in files:
            files[extension].append(i)
    return files

# private function to see if we have data
def __has_files(filelist):
    return len(filelist)

# private functions to check extension
def __is_bam(extension, filelist):
    return extension == 'BAM' and __has_files(filelist)

def __is_sam(extension, filelist):
    return extension == 'SAM' and __has_files(filelist)

def __is_sff(extension, filelist):
    return extension == 'SFF' and __has_files(filelist)

def __is_fastq(extension, filelist):
    return extension == 'FASTQ' and __has_files(filelist)

def __is_fasta(extension, filelist):
    return extension == 'FASTA' and __has_files(filelist)

def __is_qual(extension, filelist):
    return extension == 'QUAL' and  __has_files(filelist)

# public function for printing errors
def print_error(error):
    sys.exit(error)

# private function to check to see how many files we found.
#    We are limiting the numbers allowed.
def __check_size(extension,max,files):
    if len(files) > max:
        print_error("Only " + str(max) + " files of type, " + extension + ", allowed.  Found " + str(len(files)) + " instead.")
    return 1

# private function to do most of the work
def __check_files(files):
    file_object = None
    num_files = 0
    pairing = None
    if options.direction:
        pairing = options.direction

    # go through file dictionary and look for types
    # set up objects based on type
    for i in sorted(files.iterkeys()):

        if __is_bam(i,files[i]):
            __check_size(i,1,files[i])
            file_object = BamFile(files[i][0],i,options.output,pairing)
            num_files = 1
            break
        
        if __is_sam(i,files[i]):
            __check_size(i,1,files[i])
            file_object = BamFile(files[i][0],i,options.output,pairing)
            num_files = 1
            break

        if __is_sff(i,files[i]):
            __check_size(i,1,files[i])
            file_object = SffFile(files[i][0],pairing,options.output)
            num_files = 1
            break
        
        if __is_fastq(i,files[i]):
            __check_size(i,2,files[i])
            file_object = FastqFile(files[i],options.sample,options.qual_format,pairing,options.output)
            num_files = len(files[i])
            break

        # need to check if we have quals
        if __is_fasta(i,files[i]):
            __check_size(i,2,files[i])
            quals = []

            # get quals if there
            for j in sorted(files.iterkeys()):
                if __is_qual(j,files[j]):
                    quals.append(files[j])

            # check to see if we have same number of quals
            if len(quals) == len(files[i]) or not len(quals):
                file_object = FastaFile(files[i],options.sample,options.quality_score,quals,pairing, options.output)
                num_files = len(files[i])
                break
            else:
                print_error("Number of fasta files, " + str(len(files[i])) + \
                            ", does not equal number of qual files, " + str(len(quals)) + ".") 
    return file_object, num_files

def main():

    java_heap = options.java_heap

    if not re.search("G",java_heap):
        java_heap += "G"

    # get our files and objects created
    files = __initialize_files()
    files = __read_args(args,files)
    file_object, check = __check_files(files)

    # do we have a file to do something with?

    if check:
        # check our return value during conversion
        ret_value = file_object.convert_to_unmapped_bam(java_heap)
        if ret_value:
            print "Conversion successful."
            return 0
        else:
            print_error("ERROR:  File conversion failed.")
            return -1
    else:
        print_error("Can't determine file types for " + str(args))

if __name__ == "__main__":
    sys.exit(main())
