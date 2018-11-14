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

%prog [options] assembly.fasta reference.fasta

Computes assembly coverage vs a reference on the basis of shared
kmers. The kmer size used is a command line argument.

"Covered" is a measure of how many of the reference kmers are in the
assembly.  If a particular kmer occurs three times in the reference
but only twice in the assembly, that counts for 2 out of 3.  If it
occurs more times in the assembly, it doesn't count extra here. This
is total reference coverage.

"Distinct" is a measure of how many distinct kmers in the refrence are
present in the assembly (i.e., no matter how many times a given kmer
occurs in the reference, you get credit if it occurs at least once in
the assembly).  This is coverage of distinct sequence.

"OverCovered" is a count of how many kmers occur more often in the
assembly than the reference.  I.e., if a given kmer occurs twice in
the reference, but 4 times in the assembly, that counts as +2.  This
is a measure of how much in the assembly is "duplicated" from the
reference.  It does *not* include kmers in the assembly which are no
present in the reference at all, because that is....

"Novel" is a count of kmers in the assembly which are not present at
all in the reference, counting them multiple times if present multiple
times.  This is total sequence in the assembly but not in the
reference.

"NovelDistinct" is same as "Novel", but just counting each novel kmer
once...how much unique sequence is present in the assembly.

Warning: this might not play well with flattened diploid assemblies,
as variable SNPs will each contribute K bases of difference.
""")

parser.add_option('-k', default=29, type='int', help='Kmer size (default 29)')
parser.add_option('-o', help='File prefix for table output')
parser.add_option('--no_html', action="store_false", dest="html", default=True, help='Suppress html table output')

(options, args) = parser.parse_args()

if len(args) != 2:
   parser.print_help()
   exit(1)

kmer = Kmer(options.k)

print 'Kmerizing assembly'
assembly_seqs = (s.seq for s in SeqIO.parse(args[0], 'fasta'))
assembly_counts = kmer.kmer_count_seqs(assembly_seqs)

print 'Kmerizing reference'
reference_seqs = (s.seq for s in SeqIO.parse(args[1], 'fasta'))
reference_counts = kmer.kmer_count_seqs(reference_seqs)

total = 0
covered = 0
distinct = 0
distinct_covered = 0
novel = 0
novel_distinct = 0
over = 0

print 'Computing coverage'
for k, n in reference_counts.iteritems():
    distinct += 1
    total += n
    a = assembly_counts.get(k, 0)
    covered += min(a, n)
    if a:
        distinct_covered += 1
    if a > n:
        over += a - n
        


def format_row(label, n, d):
    return [label, n, d, "%.2f%%" % (n * 100.0 / d)]

headers = ['Category', 'Count', 'RefCount', 'Pct']
data = []
data.append(format_row('Covered', covered, total))
data.append(format_row('DistinctCovered', distinct_covered, distinct))
data.append(format_row('OverCovered', over, total))

for k, n in assembly_counts.iteritems():
    if k not in reference_counts:
        novel += n
        novel_distinct += 1

data.append(format_row('Novel', novel, total))
data.append(format_row('NovelDistinct', novel_distinct, distinct))

table = SimpleTable(headers, data, 'Kmer Coverage (k=' + str(options.k) + ')')

if options.o:
   table.print_output(options.o + ".kmer_coverage", options.html)
else:
   print table.to_table()


