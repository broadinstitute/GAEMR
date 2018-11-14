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
from Bio import SeqIO
from optparse import OptionParser
from gaemr.AgpFile import AgpFile
from gaemr.SimpleTable import SimpleTable
import pysam
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()

def bam_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, re.sub(' ','',value).split(','))

parser = OptionParser(usage="usage: %prog [options] contigs.fasta")

parser.add_option('--contig_table', '-c', action="store", default=None, type='string', dest='c_table',
                  help='Contig table output header (default=%default)')
parser.add_option('--scaffold_table', '-s', action="store", default=None, type='string', dest='s_table',
                  help='Scaffold table output header (default=%default)')
parser.add_option('--fragment_bam', '-f', action="callback", default=None, type='string', dest="frag_bam",
                  callback=bam_callback, help='Fragment pair bam(s) of reads aligned to input fasta.  ' + 
                  'If > 1 file, separate filenames by comma with no space (default=%default)')
parser.add_option('--jump_bam', '-j', action="callback", default=None, type='string', dest="jump_bam",
                  callback=bam_callback, help='Jump pair bam(s) of reads aligned to input fasta.  ' + 
                  'If > 1 file, separate filenames by comma with no space (default=%default)')
parser.add_option('--unpaired_bam', '-u', action="callback", default=None, type='string', dest="unpaired_bam",
                  callback=bam_callback, help='Unpaired bams (e.g. PacBio) of reads aligned to input fasta.  ' +
                  'If > 1 file, separate filenames by comma with no space (default=%default)')
parser.add_option('--taxonomy_output', '-t', action="store", default=None, type='string', dest="tax_file",
                  help='Taxonomy output from get_blast_hit_taxonomy.py (default=%default)')
parser.add_option('--agp', '-a', action="store", default=None, dest="agp",
                  help='AGP filename used to get linking information when used with contig input (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 1:
    # NEED TO HANDLE USE CASES BETTER
    parser.error("Must supply fasta file.")
    sys.exit(1)

def __get_gc_from_seq(seq):
    gc = 0
    total = len(seq)
    gc += seq.count('G') + seq.count('C')
    return (float(gc)/total)*100, gc

def __get_seq_info(fasta):
    data = {}
    gc = {}
    try:
        for f in SeqIO.parse(open(fasta,'r'),"fasta"):
            data[f.id] = {}
            data[f.id]['length'] = len(f.seq)
            data[f.id]['gc'], gc[f.id] = __get_gc_from_seq(f.seq.upper())
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1

    return data,gc

def __get_scaffolding_info(seqs,agp):
    a = AgpFile(agp)
    scaffolds = {}
    for i in seqs:
        scaffold = a.get_scaffold_from_contig_id(i)
        seqs[i]['scaffold'] = scaffold
        if scaffold not in scaffolds:
            scaffolds[scaffold] = []
        scaffolds[scaffold].append(i)
    return seqs,a,scaffolds

def __get_taxonomy(seqs, tax_file):
    taxonomy = {}
    try:
        f = open(tax_file,'r')
        for i in f.readlines():
            if re.match('^#', i):
                continue
            else:
                data = i.split(constant.TABLE_DELIMITER)
                contig = data[0]
                tax_list = data[-1].split(';')
                taxonomy[contig] = {}
                taxonomy[contig]['hit len'] = int(re.sub(",","",data[2]))
                found_genus = 0
                for j in tax_list:
                    group,name = j.split('=')
                    if group == 'genus':
                        taxonomy[contig]['genus'] = name
                        found_genus = 1
                        break
                if not found_genus:
                    taxonomy[contig]['genus'] = "NA"
        f.close()
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1
    return seqs,taxonomy

def __get_coverage_from_bam(bam_dict,s_i_dict,agp=None):
    coverage = {}

    for bam in sorted(bam_dict.iterkeys()):
        coverage[bam] = {}
        if bam_dict[bam]:
            for bam_file in bam_dict[bam]:
                samfile = pysam.Samfile(bam_file,"rb")

                for id in s_i_dict:
                    s_id = id
                    start = 1
                    stop = s_i_dict[id]['length']
                    base_coverage = 0
                    if agp:
                        s_id,start,stop = agp.get_coordinates_from_contig_id(id)
                    for pileupcolumn in samfile.pileup(s_id,start,stop):
                        base_coverage += int(pileupcolumn.n)
                    if id not in coverage[bam]:
                        coverage[bam][id] = 0
                    coverage[bam][id] += int(float(base_coverage)/((stop - start) + 1))

                samfile.close()

    return coverage

def __get_coverage_string(id,cov_dict):
    cov_string = " ("
    total = 0
    for b in sorted(cov_dict.iterkeys()):
        if cov_dict[b]:
            if id in cov_dict[b]:
                total += cov_dict[b][id]
                cov_string += str("%.0f" % cov_dict[b][id]) + "/"
            else:
                cov_string += "0/"
        else:
            cov_string += "0/"
            
    return str("%.0f" % total) + cov_string[:-1] + ")"

def __get_blast_strings(id,tax_dict,length):
    if id in tax_dict:
        return tax_dict[id]['genus'],"%.2f" % (float(tax_dict[id]['hit len'])/length * 100)
    return "NA","NA"
    
def __print_contig_stats(si_dict,cov_dict,tax_dict):
    #Contig scaffold, length, gc, cov (frg/jump) blast_hit blast_cvg ambig
    title = "Detailed Contig Stats"
    headers = ["Contig","Scaffold","Length","GC","Coverage(F/J/LR)","BLAST Hit","BLAST Covered"]
    data = []

    for c in sorted(si_dict.iterkeys()):
        tmp = []
        tmp.append(c)
        if 'scaffold' in si_dict[c]:
            tmp.append(si_dict[c]['scaffold'])
        else:
            tmp.append("NA")
        tmp.append(si_dict[c]['length'])
        tmp.append("%.2f" % si_dict[c]['gc'])
        tmp.append(__get_coverage_string(c,cov_dict))
        hit,covered = __get_blast_strings(c,tax_dict,si_dict[c]['length'])
        tmp.append(hit)
        tmp.append(covered)
        data.append(tmp)
    
    st = SimpleTable(headers,data,title)
    output = None
    if options.c_table:
        output = options.c_table + ".contig_detail"
    st.print_output(output,options.html)

def __get_scaffold_gc(s_list,si_dict,cg_dict):
    gc = 0
    total = 0
    for c in s_list:
        gc += cg_dict[c]
        total += si_dict[c]['length']
    if total:
        return "%.2f" % ((float(gc)/total) * 100)
    return 0

def __get_scaffold_coverage(s_list,cov_dict):
    cov_string = " ("
    total = 0
    total_dict = {}
    total_contigs = 0
    
    for id in s_list:
        total_contigs += 1
        for b in sorted(cov_dict.iterkeys()):
            if b not in total_dict:
                total_dict[b] = 0
            if cov_dict[b]:
                if id in cov_dict[b]:
                    total += cov_dict[b][id]
                    total_dict[b] += cov_dict[b][id]

    if total_dict:
        for i in sorted(total_dict.iterkeys()):
            cov_string += str("%.0f" % (float(total_dict[i])/total_contigs)) + "/"
    else:
        cov_string += ''.join(["0/" for i in xrange(3)])
    return str("%.0f" % (float(total)/total_contigs)) + cov_string[:-1] + ")"

def __get_longest_hit(blast):
    max = 0
    id = None
    for i in blast:
        if blast[i] > max:
            id = i
            max = blast[i]
    return id,max

def __get_scaffold_blast_hits(s_list,tax_dict,si_dict):
    blast = {}
    length = 0
    for id in s_list:
        length += si_dict[id]['length']
        if id in tax_dict:
            if tax_dict[id]['genus'] not in blast:
                blast[tax_dict[id]['genus']] = 0
            blast[tax_dict[id]['genus']] += int(tax_dict[id]['hit len'])
    genus,hit_length = __get_longest_hit(blast)
    return genus,"%.2f" % ((float(hit_length)/length) * 100)

def __print_scaffold_stats(si_dict,s_dict,cg_dict,cov_dict,tax_dict,agp_obj):
    title = "Detailed Scaffold Stats"
    headers = ["Scaffold","Num Contigs","Length","GC","Coverage(F/J/LR)","BLAST Hit","BLAST Covered","Circular"]
    data = []
    for s in sorted(s_dict.iterkeys()):
        tmp = []
        tmp.append(s)
        tmp.append(len(s_dict[s]))
        tmp.append(agp_obj.get_scaffold_length(s))
        tmp.append(__get_scaffold_gc(s_dict[s],si_dict,cg_dict))
        tmp.append(__get_scaffold_coverage(s_dict[s],cov_dict))
        hit,hit_len = __get_scaffold_blast_hits(s_dict[s],tax_dict,si_dict)
        tmp.append(hit)
        tmp.append(hit_len)
        tmp.append("NA")
        data.append(tmp)
        
    st = SimpleTable(headers,data,title)
    output = None
    if options.s_table:
        output = options.s_table + ".scaffold_detail"
    st.print_output(output,options.html)
    
    
def main():

    seq_info_dict,contig_gc_dict = __get_seq_info(args[0])
    
    agp_obj = None
    scaffold_dict = {}
    
    if options.agp:
        seq_info_dict,agp_obj,scaffold_dict = __get_scaffolding_info(seq_info_dict,options.agp)

    bam_dict = {'frag': options.frag_bam, 'jump': options.jump_bam, 'unpaired': options.unpaired_bam}
    bam_coverage_dict = __get_coverage_from_bam(bam_dict,seq_info_dict,agp_obj)

    taxonomy_dict = {}
    if options.tax_file:
        seq_info_dict,taxonomy_dict = __get_taxonomy(seq_info_dict, options.tax_file)

    
    __print_contig_stats(seq_info_dict,bam_coverage_dict,taxonomy_dict)
    if options.agp:
        __print_scaffold_stats(seq_info_dict,scaffold_dict,contig_gc_dict,bam_coverage_dict,taxonomy_dict,agp_obj)
    

    return 0
        
if __name__ == "__main__":
    sys.exit(main())



