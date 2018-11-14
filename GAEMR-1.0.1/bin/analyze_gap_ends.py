#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import sys
from itertools import *
from operator import *
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from gaemr.Kmer import Kmer
from gaemr.SimpleTable import SimpleTable
from gaemr.AgpFile import AgpFile
from gaemr.ChartUtilities import ChartUtilities
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()
import numpy as np


parser = OptionParser(usage="usage: %prog [options] <assembly.contigs.fasta> <assembly.agp>")
parser.add_option('-k', default=29, type='int', help='Kmer size (default 29)')
parser.add_option('-e', default=50, type='int', help='Gap end length (default 50)')
parser.add_option('-l', default=75, type='int', help='Percent distinct threshold for calling a gap low complexity (default 75)')
parser.add_option('-c', action="store", type='string', dest="header", default="gap_analysis", help='Header for chart output (default = %default)')
parser.add_option('-t', action="store", type='string', dest="t_header", default="gap_analysis", help='Header for table output (default = %default)')
parser.add_option('--no_html', action="store_false", dest="html", default=True, help='Suppress html table output')

(options, args) = parser.parse_args()

chart_tools = ChartUtilities()

if len(args) < 1:
    parser.error("Must supply a contigs fasta file.")
    sys.exit(-1)

def get_simple_sequences(nuc_dict,seq):
    for nuc in nuc_dict.keys():
        #print "NUC",nuc
        #print "SEQ", seq
        nuc_len = len(nuc)
        starts = [match.start() for match in re.finditer(re.compile(nuc), str(seq))]
        #print "STARTS",starts
        a = np.array(starts)
        for x in np.array_split(a,np.where(np.diff(a)!=nuc_len)[0]+1):
            if len(x) >= 20/nuc_len:
                nuc_dict[nuc] += len(x)*nuc_len
    return nuc_dict

def get_distinct_percent(count):
    distinct = 0
    total = 0
    for k,n in count.iteritems():
        distinct += 1
        total += n        
    return (float(distinct)/total)*100

def get_avg_copy_number(kmer_list,global_counts,kmer):
    sum = 0
    count = 0
    for k in kmer_list:
        rc = kmer.kmer_rc(k)
        if ((not k in global_counts) and (not rc in global_counts)):
            print "WARNING: Neither of the following kmers could be found:",kmer.kmer_to_seq(k),k,kmer.kmer_to_seq(rc),rc
            continue
        sum = sum + global_counts[k]
        count = count + 1
    if not count:
        return 0
    else:
        return float(sum)/count

def analyze_seq(seq,flank_kmer,asm_kmer,asm_counts):
    flank_counts = flank_kmer.kmer_count_seq(seq)
    gc = GC(seq)
    dist = get_distinct_percent(flank_counts)
    copy = get_avg_copy_number(asm_kmer.kmer_list(seq),asm_counts,asm_kmer)
    return (gc,dist,copy)

def __initialize_nucs():
    nucs = ["A","C","T","G","AC","AT","AG","CT", "CG", "TG", "AAC", "AAG", "AAT", "ACC", "ACG", "ACT", "AGC", "AGG", "AGT", "ATC", "ATG", "ATT", "CCG", "CCT", "CGG", "CGT", "CTG", "CTT", "GGT", "GTT"]
    cg_ss_counts = {}
    for nuc in nucs:
        cg_ss_counts[nuc] = 0;
    ue_ss_counts = {}
    for nuc in nucs:
        ue_ss_counts[nuc] = 0;

    return cg_ss_counts, ue_ss_counts

def __get_contig_seqs(file):
    # Store contig sequences
    contig_seqs = {}
    for c in SeqIO.parse(file, 'fasta'):
        contig_seqs[c.id] = c.seq
    return contig_seqs

def __get_assembly_kmer_counts(file, asm_kmer):
    # Store assembly sequences for counting
    asm_seqs = [s.seq for s in SeqIO.parse(file, 'fasta')]
    return asm_kmer.kmer_count_seqs(asm_seqs)

def __write_output_data(agp_file=None, flank_kmer=None, asm_kmer=None, contig_seqs=None, asm_counts=None, cg_ss_counts=None, ue_ss_counts=None):
    #captured gap stats
    cg_seqs = []
    cg_sizes = []
    cg_dist = []
    cg_cn = []
    cg_gc = []
    
    #uncaptured end stats
    ue_seqs = []
    ue_dist = []
    ue_cn = []
    ue_gc = []
    
    print "Pulling gap end sequence..."
    agp = AgpFile(agp_file)
    scaffolds = agp.get_agp_scaffolds()

    for scaffold in scaffolds:
        print "Scaffold:", scaffold
        for record in agp.get_agp_file_record(scaffold):
            if (not agp.is_gap(scaffold,record)):
                ctg_seq = contig_seqs[agp.get_contig_id(scaffold,record)]
                if (len(ctg_seq) < options.e):
                    print "WARNING:", agp.get_contig_id(scaffold,record), " is less than ", options.e, ". Ignoring."
                    continue
            if (record == 1):
                #print "\tFirst contig..."
                #print "SEQ: ", contig_seqs[agp.get_contig_id(scaffold,record)][:options.e]
                (gc,dist,copy) = analyze_seq(contig_seqs[agp.get_contig_id(scaffold,record)][:options.e],flank_kmer,asm_kmer,asm_counts)
                ue_dist.append(dist)
                ue_cn.append(copy)
                ue_gc.append(gc)
                ue_ss_counts = get_simple_sequences(ue_ss_counts,contig_seqs[agp.get_contig_id(scaffold,record)][:options.e])
                #print "UE", ue_ss_counts
            if (record == len(agp.get_agp_file_record(scaffold))):
                #print "\tLast contig..."
                #print "SEQ: ", contig_seqs[agp.get_contig_id(scaffold,record)][-options.e:]

                (gc,dist,copy) = analyze_seq(contig_seqs[agp.get_contig_id(scaffold,record)][-options.e:],flank_kmer,asm_kmer,asm_counts)
                ue_dist.append(dist)
                ue_cn.append(copy)
                ue_gc.append(gc)
                ue_ss_counts = get_simple_sequences(ue_ss_counts,contig_seqs[agp.get_contig_id(scaffold,record)][-options.e:])
                #print "UE", ue_ss_counts
            if (agp.is_gap(scaffold,record)):
                #print "Gap..."
                left_ctg_record = record-1
                right_ctg_record = record+1

                if (len(contig_seqs[agp.get_contig_id(scaffold,left_ctg_record)]) < options.e) or (len(contig_seqs[agp.get_contig_id(scaffold,right_ctg_record)]) < options.e):
                    print "Warning: This gap is flanked by a contig less than",options.e,"long. Skipping analysis."
                    continue
                left_seq = contig_seqs[agp.get_contig_id(scaffold,left_ctg_record)][-options.e:]
                #left_seq = contig_seqs[agp.get_contig_id(scaffold,left_ctg_record)][:options.e]
                #print "LEFT: ", left_seq
                (left_gc,left_dist,left_copy) = analyze_seq(left_seq,flank_kmer,asm_kmer,asm_counts)
                cg_ss_counts = get_simple_sequences(cg_ss_counts,left_seq)
                #print "CG", cg_ss_counts
                #print left_seq
                #print left_gc,left_dist,left_copy

                right_seq = contig_seqs[agp.get_contig_id(scaffold,right_ctg_record)][:options.e]
                #right_seq = contig_seqs[agp.get_contig_id(scaffold,right_ctg_record)][-options.e:]
                #print "RIGHT SEQ: ", right_seq
                (right_gc,right_dist,right_copy) = analyze_seq(right_seq,flank_kmer,asm_kmer,asm_counts)
                cg_ss_counts = get_simple_sequences(cg_ss_counts,right_seq)
                #print "CG", cg_ss_counts
                #print right_gc,right_dist,right_copy
                
                cg_sizes.append(agp.get_feature_length(scaffold,record))
                cg_dist.append((left_dist+right_dist)/2)
                cg_gc.append((left_gc+right_gc)/2)
                cg_cn.append((left_copy+right_copy)/2)


    headers = ["Metric", "Uncaptured Ends", "Captured Gaps"]
    table = SimpleTable(headers,[],"Gap Analysis Metrics")

    ss_table = SimpleTable(["Sequence","Uncaptured Ends Bases", "Uncaptured Ends Percent", "Captured Gap Bases", "Captured Gap Percent"],[],"Gap Simple Sequence Analysis")

    if (len(cg_sizes)>0):
        chart_tools.gen_histogram(cg_sizes,"Gap Sizes", "Number of Gaps","Histogram of Captured Gap Sizes",options.header+".cg_sizes")
        chart_tools.gen_histogram(cg_dist,"Gap Flank Complexity", "Number of Gaps","Histogram of Gap Flank Complexity",options.header+".cg_distinctness")
        chart_tools.gen_histogram(cg_cn,"Gap Flank Copy Number", "Number of Gaps","Histogram of Gap Flank Copy Number",options.header+".cg_copy_number")
        chart_tools.gen_histogram(cg_gc,"Gap Flank GC", "Number of Gaps","Histogram of Gap Flank GC",options.header+".cg_gc")

        table.add_row(["Number",len(ue_dist),len(cg_dist)])
        table.add_row(["Average Complexity","%.0f"%(sum(ue_dist)/len(ue_dist)), "%.0f"%(sum(cg_dist)/len(cg_dist))])
        table.add_row(["Less than "+str(options.l)+"% Complex",len(filter(lambda x: x < options.l,ue_dist)),len(filter(lambda x: x < options.l,cg_dist))])
        table.add_row(["Average GC", "%.0f"%(sum(ue_gc)/len(ue_gc)), "%.0f"%(sum(cg_gc)/len(cg_gc))])
        table.add_row(["Less than 30% GC", len(filter(lambda x: x < 30,ue_gc)), len(filter(lambda x: x < 30,cg_gc))])
        table.add_row(["Greater than 70% GC", len(filter(lambda x: x > 70,ue_gc)), len(filter(lambda x: x > 70,cg_gc))])
        table.add_row(["Average Copy Number", "%.0f"%(sum(ue_cn)/len(ue_cn)), "%.0f"%(sum(cg_cn)/len(cg_cn))])
        ss_table.add_row(["End Bases", len(ue_dist)*options.e, "",len(cg_dist)*options.e,""])
        ss_table.add_row(["Total SS",
                          sum(ue_ss_counts.itervalues()),
                          "%.2f%%"%(float(sum(ue_ss_counts.itervalues())*100)/(len(ue_dist)*options.e)),
                          sum(cg_ss_counts.itervalues()),
                          "%.2f%%"%(float(sum(cg_ss_counts.itervalues())*100)/(len(cg_dist)*options.e))])
        for n in cg_ss_counts:
            ss_table.add_row([n,
                          ue_ss_counts[n],
                          "%.2f%%"%(float(ue_ss_counts[n]*100)/(len(ue_dist)*options.e)),
                          cg_ss_counts[n],
                          "%.2f%%"%(float(cg_ss_counts[n]*100)/(len(cg_dist)*options.e))])
        
    else:
        print "No captured gaps."
        table.add_row(["Number",len(ue_dist),len(cg_dist)])
        table.add_row(["Average Uniqueness","%.0f"%(sum(ue_dist)/len(ue_dist)), "N/A"])
        table.add_row(["Less than "+str(options.l)+"% distinct","%.2f"%(len(filter(lambda x: x < options.l,ue_dist))),"%.2f"%(len(filter(lambda x: x < options.l,cg_dist)))])
        table.add_row(["Average GC", "%.0f"%(sum(ue_gc)/len(ue_gc)), "N/A"])
        table.add_row(["Less than 30% GC", len(filter(lambda x: x < 30,ue_gc)), len(filter(lambda x: x < 30,cg_gc))])
        table.add_row(["Greater than 70% GC", len(filter(lambda x: x > 70,ue_gc)), len(filter(lambda x: x > 70,cg_gc))])
        table.add_row(["Average Copy Number", "%.0f"%(sum(ue_cn)/len(ue_cn)), "N/A"])
        ss_table.add_row(["End Bases", len(ue_dist)*options.e, "","NA",""])
        ss_table.add_row(["Total SS",
                          sum(ue_ss_counts.itervalues()),
                          "%.2f%%"%(float(sum(ue_ss_counts.itervalues())*100)/(len(ue_dist)*options.e)),
                          "NA",
                          "NA"])
        for n in cg_ss_counts:            
            ss_table.add_row([n,
                          ue_ss_counts[n],
                          "%.2f%%"%(float(ue_ss_counts[n]*100)/(len(ue_dist)*options.e)),
                          "NA",
                          "NA"])
            
    chart_tools.gen_histogram(ue_dist,"End Distinctness", "Number of Uncaptured Ends","Histogram of End Complexity",options.header+".ue_distinctness")
    chart_tools.gen_histogram(ue_cn,"End Copy Number", "Number of Uncaptured Ends","Histogram of End Copy Number",options.header+".ue_copy_number")
    chart_tools.gen_histogram(ue_gc,"End GC", "Number of Uncaptured Ends","Histogram of End GC",options.header+".ue_gc")

    table.print_output(options.t_header + ".gap_analysis",options.html)
    ss_table.print_output(options.t_header + ".gap_ss_analysis",options.html)
    
def main():
    
    ###Analyzing Distinct 5-mers at contig ends
    print "Loading assembly kmers..."

    #Store Nucs for Simple Sequence analysis
    cg_ss_counts, ue_ss_counts = __initialize_nucs()

    # A 5-base Kmer for flanks
    flank_kmer = Kmer(5,rc=False)
    # A 29-base Kmer for assembly
    asm_kmer = Kmer(options.k)

    contig_seqs = __get_contig_seqs(args[0])

    asm_counts = __get_assembly_kmer_counts(args[0], asm_kmer)

    __write_output_data(agp_file=args[1], flank_kmer=flank_kmer, asm_kmer=asm_kmer, contig_seqs=contig_seqs, \
                        asm_counts=asm_counts,cg_ss_counts=cg_ss_counts, ue_ss_counts=ue_ss_counts)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
