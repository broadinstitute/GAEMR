#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import os
import shutil
import sys
import re
from multiprocessing import Process,Queue,Manager
import multiprocessing
from optparse import OptionParser
import gaemr.PlatformConstant as pc
from gaemr.WrapperUtils import *

constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option('--reference', '-r', action="store", type='string',default=None,dest="reference",
    help='Reference file to use for analysis (default=%default)')
parser.add_option('--scaffolds', '-s', action="store", type='string',default=None,dest="scaffolds",
    help='Scaffold fasta file to use for analysis (default=%default)')
parser.add_option('--contig', '-c', action="store", type='string',default=None,dest="contigs",
    help='Contig fasta file to use for analysis (default=%default)')
parser.add_option('--agp', '-a', action="store", type='string',default=None,dest="agp",
    help='Assembly AGP file to use for analysis (default=%default)')
parser.add_option('--read_list','-l', action="store", type='string',default=None,dest="read_list_file",
    help='File containing the list of reads to use for alignment (default=%default)')
parser.add_option('--threads','-t', action="store", type='int',default=1,dest="threads",
    help='Number of threads to use for multi-threadable processes (default=%default)')
parser.add_option('--window_size','-w', action="store", type='int',default=1000,dest="window_size",
    help='Window size for various chart outputs (default=%default)')
parser.add_option('--assembly_name','-n',action="store",type='string',default="assembly",dest="assembly_name",
    help='Name of assembly (default=%default)')
parser.add_option('--assembler','-m',action="store",type='string',default="assembler",dest="assembler",
    help='Assembler used to generate consensus (default=%default)')
parser.add_option('--output','-o',action="store",type='string',default="assembly",dest="output",
    help='Output header (default=%default)')
parser.add_option('--analyze_rna','-g',action="store_true",default=None,dest='analyze_rna',
    help='Run rna analysis on assembly (default=%default)')
parser.add_option('--blast_task','-b',action="store",default='megablast',dest='blast_task',
    help='Type of blast to run for contamination check \n blastn, dc-megablast, megablast (default=%default)')
parser.add_option('--aligner','-i',action="store",default='bwa',dest='aligner',
    help='Aligner to use for reads alignments (bwa or bowtie2) (default=%default)')
parser.add_option('--kmer_size','-k',action="store",type=int,default=29,dest="kmer_size",
    help='k-mer size for repeat and k-mer coverage stats (default=%default)')
parser.add_option('--no_blast',action="store_true",default=False,dest='no_blast',
    help='Do not run blast jobs (default=%default)')
parser.add_option('--force',action="store_true",default=False,dest='force',
    help='Force remove a previous gaemr dir (default=%default)')
parser.add_option('--minScaffSize', action="store", type=int, default=1, dest='minScaffSize',
    help='Minimum scaffolds size for inclusion in analysis (default=%default)')
parser.add_option('--minConSize', action="store", type=int, default=1, dest='minConSize',
    help='Minimum contig size for inclusion in analysis (default=%default)')
parser.add_option('--minGapSize', action="store", type=int, default=10, dest='minGapSize',
    help='Minimum length of run of Ns to determine a gap (default=%default)')
parser.add_option('--base_url', action="store", type='string', default=None, dest='base_url',
    help='Home url for linking back from GAEMR website (default=%default)')

(options, args) = parser.parse_args()

# check options and args; determine analysis types based on that.
# make an analysis dir (GAEMR)
# make a work, table and chart dir.  work will get intermediate files, table will get tables, and charts will get charts.
# Major Analysis:
#  1. Basic Assembly Stats w/ Cumulative Graphs - done
#  2. Run BLAST and taxonomy - done
#  2a.  RNA - done
#  3. Align reads to assembly
#     a.  Plot coverage - need to plot window coverage
#     b.  Plot insert sizes - done; need aligned bams
#     c.  Get general alignment stats - done; need aligned bams
#  4. K-mer Analysis
#  5. If ref, reference analysis stuff
#     a. nucmer
#     b. % covered
#     c. align reads?
#     d. scaffold accuracy - Terry
#  6. Make detailed table - needs to happen at end

if not options.contigs and not options.scaffolds:
    parser.error("Must supply scaffold fasta using -s, or a contig fasta file using -c (agp optional).")
    sys.exit(1)

def run_blast_analysis(wu=None, query=None, threads=None, agp=None, blast_task=None):
    blast_xml=wu.run_blast(query=query, threads=threads, blast_task=blast_task)
    parsed_blast=wu.parse_blast_xml(blast_xml=blast_xml)
    taxonomy_heatmap=wu.get_blast_hit_taxonomy(parsed_blast=parsed_blast, query=query)
    wu.blast_map(taxonomy_heatmap=taxonomy_heatmap, agp=agp)
    #wu.blast_map(taxonomy_heatmap(wu.get_blast_hit_taxonomy(parsed_blast=wu.parse_blast_xml(blast_xml=wu.run_blast(query=query, threads=threads)), query=query)), agp=agp)

def nucmer_analysis(wu=None, ref=None, query=None):
    coords_file = wu.run_nucmer(ref=ref, query=query)
    wu.compare_to_reference(coords_file=coords_file)

def __align_reads(wu=None, ar_dict=None, read_file_dict=None, db=None, db_header=None):

    #{'All.PacBio': {'unpaired': {'files': '/gsap/assembly_analysis/sean/GAEMR_TEST_DIR/gaemr/work/All.PacBio.unmapped.bam', 'length': '900', 'insert_size': ''}},
    # 'C0DJUACXX.2.Solexa-76388': {'jump': {'files': '/gsap/assembly_analysis/sean/GAEMR_TEST_DIR/gaemr/work/C0DJUACXX.2.Solexa-76388.unmapped.bam', 'length': '101', 'insert_size': '3000'}},
    # 'C0DJUACXX.1.Pond-121053': {'fragment': {'files': '/gsap/assembly_analysis/sean/GAEMR_TEST_DIR/gaemr/work/C0DJUACXX.1.Pond-121053.unmapped.bam', 'length': '101', 'insert_size': '200'}}}
    aligner=options.aligner
    threads=options.threads

    make_index = 1
    for group in read_file_dict:
        for type in read_file_dict[group]:
            if type not in ar_dict:
                ar_dict[type] = ''
            file = read_file_dict[group][type]['file']
            output_header = re.sub("unmapped.bam",db_header,file)
            i_size = read_file_dict[group][type]['insert_size']
            align_string = '-s'
            if int(read_file_dict[group][type]['length']) > 200:
                align_string = '-l'
                if type != 'unpaired':
                    aligner = 'BowTie'
            wu.align_reads(unmapped_bam=file, ref=db, threads=threads, aligner=aligner, ref_header=db_header, make_index=make_index, align_type=align_string)
            make_index = 0
            ar_dict[type] += 'group=' + group + ';file=' + output_header + '.bam;insert_size=' + i_size + ";"

def __get_aligned_read_files(a_dict):
    files = {}

    if a_dict:
        for type in a_dict.keys():
            files[type] = {}
            tmp = a_dict[type].split(';')
            group = ''
            for i in tmp[:-1]:
                value = i.split('=')
                if re.match('group',value[0]):
                    group = value[1]
                    files[type][group] = {}
                    continue
                files[type][group][str(value[0])] = value[1]
    return files

def __analyze_read_alignment_data(wu=None, aligned_reads_dict=None, db=None, db_header=None, output=None):
    files = __get_aligned_read_files(aligned_reads_dict)
    bam_files = []

    for type in files.keys():
        is_files = []

        for group in files[type].keys():

            bam_files.append(files[type][group]['file'])
            if type == 'unpaired':
                continue
            ins_size_file = wu.run_insert_size(bam_file=files[type][group]['file'], insert_size=files[type][group]['insert_size'])
            is_files.append(ins_size_file)

        if is_files:
            wu.plot_insert_size(insert_size_files=is_files, output_base=type, ref_header=db_header)

    wu.get_simple_bam_stats(bam_files=bam_files, name=output, ref_header=db_header)
    wu.get_bam_coverage_stats(bam_files=bam_files, name=output, ref_header=db_header)
    wu.get_bam_coverage_stats(bam_files=bam_files, name=output, ref_header=db_header, want_phys_cvg='True')
    wu.generate_bam_plots(bam_files=bam_files, ref=db, name=output, ref_header=db_header, window_size=options.window_size)
    wu.identify_coverage_anomalies(bam_files=bam_files, name=output, ref_header=db_header, window_size=options.window_size)

def __make_detailed_table(wu=None, aligned_reads_dict=None, contigs=None, agp=None, taxonomy_output=None):
    files = None
    if aligned_reads_dict:
        files=__get_aligned_read_files(aligned_reads_dict)
    
    contig_detail = wu.make_detailed_table(aligned_bam_dict=files, contigs=contigs, agp=agp, taxonomy_output=taxonomy_output)
    wu.blast_bubbles(taxonomy_output=taxonomy_output, contig_detail=contig_detail)

def __run_analysis():
    manager = Manager()

    wu = WrapperUtils(force=options.force, output_header=options.output, threads=options.threads)

    scaffolds, contigs, agp = wu.standardize_file_inputs(scaffolds=options.scaffolds, contigs=options.contigs,
                                                         agp=options.agp, minScaffSize=options.minScaffSize,
                                                         minConSize=options.minConSize, minGapSize=options.minGapSize)

    [ref,read_list_file] = wu.copy_files(wu.get_work_dir(),[options.reference,options.read_list_file])

    os.chdir(wu.get_gaemr_dir())


    bas = Process(target=wu.get_basic_assembly_stats,args=(),kwargs={'assembler':options.assembler,
                                                                     'contigs':contigs, 'agp':agp})
    bas.start()

    age = Process(target=wu.analyze_gap_ends,args=(),kwargs={'contigs':contigs, 'agp':agp})
    age.start()

    rb = ar = aar = None

    if not options.no_blast:
        rb = Process(target=run_blast_analysis, args=(),kwargs={'wu':wu, 'query':contigs,'threads':options.threads, 'agp':agp, 'blast_task':options.blast_task})
        rb.start()

    if options.analyze_rna:
        aar = Process(target=wu.analyze_assembly_rna,args=(),kwargs={'query':contigs})
        aar.start()

    bas.join()
    age.join()

    if aar:
        aar.join()
    if rb:
        rb.join()

    kcn = Process(target=wu.get_kmer_copy_number,args=(), kwargs={'fasta':scaffolds, 'kmer_size':options.kmer_size})
    kcn.start()
    kcn.join()

    ar = arad = None

    aligned_reads_dict = None
    
    if options.read_list_file:
        # get read files for alignment
        read_file_dict = wu.get_read_files(read_list_file=read_list_file)
        reads_to_align = wu.check_read_files(read_file_dict=read_file_dict)
        
        aligned_reads_dict = manager.dict()
        ar = Process(target=__align_reads, args=(), kwargs={'wu':wu, 'ar_dict':aligned_reads_dict, 'read_file_dict':reads_to_align,
                                                            'db':scaffolds, 'db_header':'scaffolds'})
        ar.start()
        ar.join()

        arad = Process(target=__analyze_read_alignment_data,args=(), kwargs={'wu':wu, 'aligned_reads_dict':aligned_reads_dict,
                                                                             'db':scaffolds, 'db_header':'scaffolds', 'output':options.output})
        arad.start()
        arad.join()

        

    if ref and os.path.exists(ref):
        rn = Process(target=nucmer_analysis,args=(), kwargs={'wu':wu, 'ref':ref, 'query':scaffolds})
        rn.start()
        rn.join()
        
        kmc = Process(target=wu.run_kmer_coverage,args=(), kwargs={'ref':ref, 'query':scaffolds, 'kmer_size':options.kmer_size})
        kmc.start()
        kmc.join()
        
        sa = Process(target=wu.run_scaffold_accuracy,args=(), kwargs={'ref':ref, 'query':scaffolds})
        sa.start()
        sa.join()

        ar = arad = None
        if options.read_list_file:
            ref_aligned_reads_dict = manager.dict()
            ar = Process(target=__align_reads, args=(), kwargs={'wu':wu, 'ar_dict':ref_aligned_reads_dict, 'read_file_dict':reads_to_align,
                                                                'db':ref, 'db_header':'reference'})
            ar.start()
        if ar:
            ar.join()
            arad = Process(target=__analyze_read_alignment_data,args=(), kwargs={'wu':wu, 'aligned_reads_dict':ref_aligned_reads_dict,
                                                                                 'db':ref, 'db_header':'reference', 'output':options.output})
            arad.start()
            arad.join()

    taxonomy_output = os.path.join(wu.get_table_dir(), wu.get_output_header() + '.blast_hit_taxonomy.table.txt')
    if not os.path.exists(taxonomy_output):
        taxonomy_output=None

    mdt = Process(target=__make_detailed_table,args=(), kwargs={'wu':wu, 'aligned_reads_dict':aligned_reads_dict, 'contigs':contigs,
                                                                'agp':agp, 'taxonomy_output':taxonomy_output})
    mdt.start()
    mdt.join()

    wu.generate_html(base_url=options.base_url)
    
    return wu.check_pids()

def main():

    return __run_analysis()

if __name__ == "__main__":
    sys.exit(main())
