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
import random
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options] <sequence.fasta>")
parser.add_option('-r', default=0, type='int', help='Randomly sample this many sequences.')
parser.add_option('-l', default=0, type='int', help='Subset this many of the largest sequences')
#parser.add_option('-s', action="store", type='string', dest="seqs_file", default=None, help='Subset sequences with specified list of headers.')
parser.add_option('-o', action="store", type='string', dest="output", default="subset.fasta", help='Name of output file.')
(options, args) = parser.parse_args()

if len(args) < 1:
    parser.error("Must supply a contigs fasta file.")
    sys.exit(-1)
    

fasta = args[0]
sequences = []
count = 0
outFile = open(options.output,"w")

print "Opening", fasta
for index, record in enumerate(SeqIO.parse(fasta,"fasta")):   
    id = record.id
    length = len(record.seq)
    sequence = record.seq
    
    if (index and index%10000 == 0):
        print "Parsed ", index, "sequences."
    
    sequences.append([id,length,sequence])


if (options.l):
    count = options.l
    print "Sorting sequences to select longest", count, "..."
    for id, length, seq in ( sorted (sequences, key=lambda sequence: sequence[1], reverse=True)):
        print id,length
        if (count):
            outFile.write(">"+id+"\n")
            outFile.write(str(seq)+"\n")
            count = count - 1
        else:
            break

if (options.r):
    count = options.r
    print "Selecting", count, "random sequences..."
    while count:
        rand = random.randint(0,len(sequences)-1)
        (id,length,seq) = sequences[rand]
        del sequences[rand]
        #print "Random sequence", count, "is", rand,"!", len(sequences) , "left to select from."
        outFile.write(">"+id+"\n")
        outFile.write(str(seq)+"\n")
        count = count - 1
        
#if (options.s):
#    print "Extracted specifed sequences in", options.s, "..."
    
    

    

