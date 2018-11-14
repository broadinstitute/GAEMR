#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

#TO DO: -g is an option, but for curent design its needed. So it should be made into
# a requirement.

#usage:
#python generate_bam_plots.py ../testdata/gc_cvrg_tool/jump_reads_orig.scaffolds.sorted.bam -g ../testdata/gc_cvrg_tool/submission.assembly.fasta -o gradtest

import sys
import argparse
from gaemr.BamCoverage import BamCoverage
from gaemr.AnalyzeGc import AnalyzeGc
from gaemr.SimpleFastaFile import SimpleFastaFile
from gaemr.ChartUtilities import *
from gaemr.Chart import *
import re


#parse arguments
parser = argparse.ArgumentParser(description="Generate a set of plots for coverage and GC% over a set of sliding windows")
parser.add_argument("--gc", "-g",type=str,help="please provide a reference file to calculate GC percentage")
parser.add_argument("--dump_data", "-d",type=str,help="dump data structure to file: please provide a file name for this option")
parser.add_argument("bamfiles",help="at least one BAM file, or a list of space separated BAM files", nargs="+")
parser.add_argument('--output', '-o',type=str,help="output header name")
parser.add_argument('--window_size','-w', default=1000,type=int,help='window size for coverage analysis')
parser.add_argument('--histograms','-hi', type=str, help='histogram output name (optional)')
parser.add_argument('--lines','-l', type=str, help='line plot of all BAMs and GC plot (multi axis plot for scale)')
parser.add_argument('--heatgc','-hgc', type=str, help='coverage and GC relationship as a line and heat map plot')
parser.add_argument('--selected','-s', type=int, help='show only this many (n) scaffolds. Note: this will show the first n scaffolds of each bam')
args = parser.parse_args()


#process user entered BAM files to generate coverage plots. Return a list of
#BamCoverage objects
def process_bams_for_coverage():
    """iterate through arg entered BAMs"""
    processed_bams=[]
    for bam in args.bamfiles:
        print "processing bam " + bam + " for coverage .... "
        processed_bams.append(BamCoverage(bam,
                                          args.output,
                                          args.window_size
                                          ))
    return processed_bams
    



#get the file name out of a full dirctory path
def get_file_name_out_of_path(name_with_path):
    """get the file name out of the path"""
    fields = name_with_path.split("/")
    return fields[-1]


#copy over coverage information to internal data structure
def transcribe_processed_cvrg(coverage_processed_bams):
    """add processed coverage to unified record"""
    coverage_by_bam=dict()
    
    for bam in coverage_processed_bams:
        coverage_by_fasta_inside_bam=dict()
        cov_table = bam.get_window_coverage_table()
        for fasta_name in sorted(cov_table.iterkeys()):
            coverage_by_fasta_inside_bam[fasta_name]=[]
            seq_windows = cov_table[fasta_name]
            for seq_window in seq_windows:
                coverage_by_fasta_inside_bam[fasta_name].append(float(seq_window.value)/float(args.window_size))
        #coverage_by_bam[get_file_name_out_of_path(bam.output_prefix)] = coverage_by_fasta_inside_bam
        coverage_by_bam[get_file_name_out_of_path(bam.bam)] = coverage_by_fasta_inside_bam 

    return coverage_by_bam


#process reference for GC% analysis
def process_ref_for_gc():
    """process for GC%"""
    processed_gc = AnalyzeGc(SimpleFastaFile(args.gc),
                             "output_gc",
                             int(args.window_size))
    return processed_gc.window_gc



#consolidate data structure(s)
def __consolidate_data__(coverage_by_bam,processed_gcs):
    consolidated_on_position = []
    headers = ['ID','Window_Number', 'Window_Start', '%GC']
    for name in sorted(processed_gcs.iterkeys()):
        gcs = processed_gcs[name]
        for window_id in sorted(gcs.iterkeys()):
            gc = gcs[window_id]
            consolidated_on_position.append([name, str(window_id), str(window_id * args.window_size), gc])

    for bam_name,fastas in coverage_by_bam.iteritems():
        headers.append(bam_name)
        window_id=0
        for fasta_id in sorted(fastas.iterkeys()):
            cvrg_values = fastas[fasta_id]
            for cvrg_value in cvrg_values:
                consolidated_on_position[window_id] += [cvrg_value]
                window_id += 1
    __dump_data_structure__(headers,consolidated_on_position)



#print data structure out to file
def __dump_data_structure__(headers, consolidated_on_position):
    out = open(args.dump_data,"w")
    out_string = ''.join([str(x) + "\t" for x in headers])
    out.write(out_string[:-1] + '\n')
    for i in consolidated_on_position:
        out_string = ''.join([str(x) + "\t" for x in i])
        out.write(out_string[:-1] + '\n')
    out.close()


    
#app starts here, need to break up main into smaller parts. too big right now.
def main():
    coverage_processed_bams = process_bams_for_coverage()
    coverage_by_bam = transcribe_processed_cvrg(coverage_processed_bams)

    processed_gc=dict()
    if args.gc:
        processed_gc = process_ref_for_gc()

    chart_tools = ChartUtilities()
    all_chart_data=[]

    #add gc data: heat map requires a 2D arrary, so tricking single [] into a
    #2D arrary by making a list of lists. internal list is only one number,1 single GC%
    gc_data=[]
    for fasta_name in sorted(processed_gc.iterkeys()):
        gcs = processed_gc[fasta_name]
        converted_to_int = []
        for window_id in sorted(gcs.iterkeys()):
            gc = gcs[window_id]
            temp=[]
            temp.append(int(gc))
            converted_to_int.append(temp)
        #print "adding " + fasta_name + " to gc data set"
        for g in converted_to_int:
            gc_data.append(g)

    all_chart_data.append(gc_data)
    
    #add cvrg data
    all_chart_data.append(coverage_by_bam)
    
    chrt = chart_tools.gen_line_and_heatmap_in_same_chart_integrated(all_chart_data,' Coverage in sequence windows of '+str(args.window_size))
    chrt.save_as(args.output + ".gc_cvg.png")

    if args.dump_data:
        __consolidate_data__(coverage_by_bam,processed_gc)

    if args.histograms:        
        obj_hist = chart_tools.gen_histograms_in_same_chart(coverage_by_bam)
        obj_hist.save_as(args.histograms + ".histogram.png")

    return 0


def add_commas_to_nums(num):
    formatted_string = ''
    data = num.split(' ')
    for i in xrange(len(data)):
        value = data[i]
        if re.match('\d+', str(value)) and not re.search('\D+', str(value)):
            if re.search('\.', str(value)):
                value = "{:,.2f}".format(float(value))
            else:
                value = "{:,}".format(int(value))
        formatted_string += value + ' '
    return formatted_string[:-1]




#app starts here, need to break up main into smaller parts. too big right now.
def main_alt():
    coverage_processed_bams = process_bams_for_coverage()
    coverage_by_bam = transcribe_processed_cvrg(coverage_processed_bams)

    processed_gc=dict()
    if args.gc:
        processed_gc = process_ref_for_gc()

    chart_tools = ChartUtilities()

    #subset of fastas (if not total) to chart
    num_fastas_to_use=0
    if args.selected:
        num_fastas_to_use = args.selected 
    else:
        num_fastas_to_use = len(processed_gc.keys())

    #collect the subset to chart
    #add gc data: heat map requires a 2D arrary, so tricking single [] into a
    #2D arrary by making a list of lists. internal list is only one number,1 single GC%
    gc_data=[]
    fasta_count=0
    for fasta_name in sorted(processed_gc.iterkeys()):
        if fasta_count<num_fastas_to_use:
            gcs = processed_gc[fasta_name]
            converted_to_int = []
            for window_id in sorted(gcs.iterkeys()):
                gc = gcs[window_id]
                temp=[]
                temp.append(gc)
                converted_to_int.append(temp)

            print "adding " + fasta_name + " to gc data set"
            for g in converted_to_int:
                gc_data.append(g)
            fasta_count+=1
        else:
            break

    if args.heatgc:
        for bam_name,bam_cvrg_data in coverage_by_bam.iteritems():
            chrt = chart_tools.gen_line_and_heatmap_in_same_chart_integrated_single(gc_data, {bam_name:bam_cvrg_data},' Coverage in sequence windows of '+str(args.window_size))
            chrt.save_as(args.output + '.' + bam_name+ ".gc_cvg.png")

    if args.dump_data:
        __consolidate_data__(coverage_by_bam,processed_gc)

    if args.histograms:        
        obj_hist = chart_tools.gen_histograms_in_same_chart(coverage_by_bam)
        obj_hist.save_as(args.histograms + ".histogram.png")



    #list of bam names
    bam_names=[]
    #coordinates where scaffold (or contig) breaks happens
    scaffold_break_coords = []
    #scaffold (or contig) names to label break lines
    scaffold_break_names = []
    #minor X tick names
    minor_xtick_labels=[]
    #minor X tick coords
    minor_xtick_coords=[]
    minor_xtick_interval=100000

    bam_index=0
    #list of lists of BAM coverages. Each index is a BAM. In each index is a list of BAM coverages
    bam_coverage_values=[]
    for bam,coverage_by_fasta in coverage_by_bam.iteritems():
        pos=0
        bam_names.append(bam)
        bam_coverage_values.append([])
        fasta_ids = coverage_by_fasta.keys()
        fasta_ids.sort()
        fasta_count=0
        for fasta_id in fasta_ids:
            if fasta_count < num_fastas_to_use:
                coverages = coverage_by_fasta[fasta_id]
                scaffold_break_names.append(fasta_id)
                scaffold_break_coords.append(pos)
                xtick_pos=0
                for cvrg in coverages:
                    bam_coverage_values[bam_index].append(cvrg)
                    if (xtick_pos*args.window_size) % minor_xtick_interval ==0:
                        if xtick_pos != 0:
                            minor_xtick_coords.append(pos)
                            minor_xtick_labels.append(add_commas_to_nums(str(xtick_pos * args.window_size)))
                    pos+=1
                    xtick_pos+=1
                fasta_count+=1
            else:
                break
        bam_index+=1

        #thin out the minor tick labels so they show up only one in some n-amount of ticks
        label_interval=5
        for i in range(len(minor_xtick_labels)):
            if i % label_interval != 0:
                minor_xtick_labels[i]=' '
            
        if args.lines:
            obj_lines  = chart_tools.gen_line_plot_with_two_y_scales(bam_coverage_values,bam_names,"","Coverage",gc_data,'GC%',"GC","Coverage and GC% Relationship","lines.gc_cvg",scaffold_break_coords, scaffold_break_names,True)
            obj_lines.save_as(args.lines + ".lines.gc_cvg.png")



        #default plot
        obj_lines = chart_tools.gen_line_plot_with_layout_managed_gc_plot(bam_coverage_values,bam_names,"Window","Coverage",gc_data,'GC%',"GC","Coverage and GC% Relationship","lines.gc_cvg",scaffold_break_coords, scaffold_break_names,minor_xtick_coords,minor_xtick_labels,True)
        obj_lines.save_as(args.output + ".gc_cvg.png")
            



    return 0



if __name__ == "__main__":
    sys.exit(main_alt())
