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
import argparse
from gaemr.BamCoverage import BamCoverage
from gaemr.SimpleStat import *
from gaemr.ChartUtilities import ChartUtilities
from gaemr.Chart import Chart
import numpy as np

#parse arguments
parser = argparse.ArgumentParser(description="a detailed report on coverage and predict regions of interest for closer analysis. Will provide a scatter plot of position to coverage of regions where the coverage is above or below three standard deviations of mean")
parser.add_argument("bamfiles",help="one BAM file",nargs="+")
parser.add_argument("output",help="a name for output files")
parser.add_argument('--window_size','-w',default=40, type=int, help='window size for coverage analysis DEFAULT=40')
parser.add_argument('--cluster_distance','-d',default=100, type=int, help='allowed distance between anomalies for clustering DEFAULT=100')
parser.add_argument('--silent','-s',action='store_true',help='do not print out pogress messages (silent run)')
args = parser.parse_args()


class IdentifyCoverageAnomalies():
    """To find anormalies in coverage and predict problems"""
    #initialize
    def __init__(self, cvrg_by_bam):
        """class initialization"""
        self.coverage_by_bam = cvrg_by_bam


    def find_regions_below_or_above_mean(self):
        coords=[]
        for bam_name,fastas in self.coverage_by_bam.iteritems():
            for fasta_name,fasta_vals in fastas.iteritems():
                stat = SimpleStat(fasta_vals)
                mean = float(stat.get_mean())
                std = float(stat.get_standard_deviation())
                pos=1
                for val in fasta_vals:
                    if float(val) > mean+(3*std) or float(val) < mean-(3*std):
                        coord=[]
                        coord.append(pos)
                        coord.append(val)
                        coords.append(coord)
                    pos+=1
        return coords



#find walls of coverages dips
#x=pos, y=cvrg,z=deltawithlastcvrg
def coverage_delta_plots(coverage_by_bam,output_file_name):
    prev=0
    x=[]
    z=[]
    pos=0
    for bam_name,fastas in coverage_by_bam.iteritems():
        for fasta_name,fasta_vals in fastas.iteritems():
            for val in fasta_vals:
                x.append(pos)
                z.append(int(val)-int(prev))
                pos+=1
        prev=0
        pos=0
    cu = ChartUtilities()
    cu. gen_scatter_plot(x,
                         z,
                         "Position",
                         "Coverage Change",
                         get_file_name_out_of_path(args.output)+" Change In Coverage",
                         output_file_name + ".coverage_change_across_genome_scatter",
                         [],
                         [])
    cu. gen_histogram(z,
                      "Change In Coverage Between Two Sequence Windows",
                      "Frequency",
                      get_file_name_out_of_path(args.output)+" Coverage Difference Between Windows Histogram",
                      output_file_name + ".coverage_change_across_genome_histogram")



#get the file name out of a full dirctory path
def get_file_name_out_of_path(name_with_path):
    """get the file name out of the path"""
    fields = name_with_path.split("/")
    return fields[-1]


#process user entered BAM files to generate coverage plots. Return a list of
#BamCoverage objects
def process_bams_for_coverage():
    """iterate through arg entered BAMs"""
    processed_bams=[]
    for bam in args.bamfiles:
        processed_bams.append(BamCoverage(bam,
                                          args.output,
                                          args.window_size
                                          ))
    return processed_bams
    

#copy over coverage information to internal data structure
def transcribe_processed_cvrg(coverage_processed_bams):
    """add processed coverage to unified record"""
    coverage_by_bam=dict()
    coverage_by_fasta_inside_bam=dict()
    
    for bam in coverage_processed_bams:
        for fasta_name,seq_windows in bam.get_window_coverage_table().iteritems():
            coverage_by_fasta_inside_bam[fasta_name]=[]
            for seq_window in seq_windows:
                coverage_by_fasta_inside_bam[fasta_name].append(float(seq_window.value)/float(args.window_size))
        coverage_by_bam[get_file_name_out_of_path(bam.output_prefix)] = coverage_by_fasta_inside_bam   
        
    return coverage_by_bam


#generate a scatter plot of position to problem in coverage
def gen_scatter_plot(coords):
    """genrate scatter plot"""
    x=[]
    y=[]
    for coord in coords:
        x.append(coord[0])
        y.append(coord[1])
    chart_tools=ChartUtilities()
    chrt = chart_tools.gen_scatter_plot(
                                        x,
                                        y,
                                        "Window Position",
                                        "Coverage",
                                        get_file_name_out_of_path(args.output)+" Coverage Above/Below 3xStDevs",
                                        args.output+"_cvrg_anomalie_scatter",
                                        [],
                                        [])
    
#to cluster windows and identify regions of problamatic coverage
def cluster(coords):
    last_coord=coords[0][0]
    clusters=[]
    curr_cluster=[]
    for coord in coords:
        if coord[0] - int(last_coord) < args.cluster_distance:
            curr_cluster.append(coord)
        else:
            clusters.append(curr_cluster)
            last_coord = coord[0]
            curr_cluster=[]
    return clusters

#generate a bar plot showing length of clusters of anormalous coverage clusters across the genome
def bar_plot_cluster_lengths(clusters):
    clus_positions=[]
    lengths=[]
    for clus in clusters:
        if len(clus)>0:
            clus_positions.append(clus[0][0])
            lengths.append(len(clus))
            
    chart_tools=ChartUtilities()
    name = chart_tools.gen_simple_bar_chart(
                                            np.arange(len(lengths)),
                                            lengths,
                                            "Cluster Group",
                                            "Cluster Size",
                                            get_file_name_out_of_path(args.output)+" Anormalous Region Cluster Length",
                                            args.output+"_cvrg_anomalie_cluster_length_bar")
    return name
        
    
    

#app starts here, need to break up main into smaller parts. too big right now.
def main():
    if not args.silent:
        print "processing bams for coverage ...."
    coverage_processed_bams = process_bams_for_coverage()
    coverage_by_bam = transcribe_processed_cvrg(coverage_processed_bams)
    analyze_tool = IdentifyCoverageAnomalies(coverage_by_bam)
    if not args.silent:
        print "finding anomalous regions of coverage ...."
    coords= analyze_tool.find_regions_below_or_above_mean()
    if not args.silent:
        print "plotting ...."
    scatter_plot_name = gen_scatter_plot(coords)
    if not args.silent:
        print "clustering ...."
    clusters = cluster(coords)
    if not args.silent:
        print "plotting cluster lengths ...."
    chart_name = bar_plot_cluster_lengths(clusters)
    if not args.silent:
        print "plotting coverage delta plot ...."
    coverage_delta_plots(coverage_by_bam,args.output)

if __name__ == "__main__":
    sys.exit(main())
