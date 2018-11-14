#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

'''
Created on Apr 30, 2012

@author: bruce
'''
       
import sys
from optparse import OptionParser
from Bio import SeqIO
from gaemr.Kmer import Kmer
from gaemr.SimpleTable import SimpleTable

parser = OptionParser(usage="""                                                                                         
usage: %prog [options] input.fasta                                                                 
                                                                                                                        
Reports copy number of kmers in the input fasta file                                                                                                                        
""")

parser.add_option('-k', default=29, type='int', help='Kmer size (default 29)')
parser.add_option('-o', help='File prefix for table output')
parser.add_option('--no_html', action="store_false", dest="html", default=True, help='Suppress html table output')

(options, args) = parser.parse_args()

if len(args) != 1:
   parser.print_help()
   exit(1)

kmer = Kmer(options.k)

seqs = (s.seq for s in SeqIO.parse(args[0], 'fasta'))
counts = kmer.kmer_count_seqs(seqs)

total = 0
distinct = 0
by_count = {}
for n in counts.itervalues():
    by_count[n] = by_count.get(n, 0) + 1
    distinct += 1
    total += n

cn = by_count.keys()
cn.sort()

mers = "%d-mers" % options.k
print '#', total, 'total', mers
print '#', distinct, 'distinct', mers

table = SimpleTable(['CopyNumber', 'Count', 'BasePct', 'BaseTotal', 'KmerPct', 'KmerTotal'], [], 'Kmer Copy Number (k=' + str(options.k) + ')')

btotal = 0.0
ktotal = 0.0
for n in cn:
    c = by_count[n]
    bpct = n * c * 100.0 / total
    btotal += bpct
    kpct = c * 100.0 / distinct
    ktotal += kpct
    table.add_row([n, c, "%.2f%%" % bpct, "%.2f%%" % btotal, "%.2f%%" % kpct, "%.2f%%" % ktotal])

if options.o:
   table.print_output(options.o + ".kmer_copy_number",options.html)
else:
   print table.to_table()

