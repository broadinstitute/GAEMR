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
from gaemr.Assembly import Assembly
from gaemr.AssemblyStatsUtil import AssemblyStatsUtil
from gaemr.SimpleTable import SimpleTable
import gaemr.PlatformConstant as pc

constant = pc.PlatformConstant()
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

parser = OptionParser(usage="usage: %prog [options] scaffolds.fasta")

parser.add_option('--assembly_name', '-n', action="store", default=None, type='string',dest="name",
                  help='An assembly id to associate with assembly sequence (default=%default)')
parser.add_option('--assembler', '-a', action="store", default=None, type='string',dest="assembler",
                  help='An assembler name to associate with assembly sequence (default=%default)')
parser.add_option('--min_contig', '-c', action="store", default=constant.MIN_CONTIG_SIZE, type='int', dest="minContig",
                  help='Minimum contig length to include in stats (default=%default)')
parser.add_option('--min_scaffold', '-s', action="store", default=constant.MIN_SCAFFOLD_SIZE, type='int', dest="minScaffold",
                  help='Minimum scaffold length to include in stats (default=%default)')
parser.add_option('--min_gap_size', '-g', action="store", default=constant.MIN_GAP_SIZE, type='int', dest="minGap",
                  help='Minimum length of stretch of Ns to consider gap sequence (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
                  help='Output table file header name (default=%default)')
parser.add_option('--cumulative_contigs', '-C', action="store_true", default=None,dest="cumulative_c",
                  help='Generate cumulative contigs output (default=%default)')
parser.add_option('--cumulative_scaffolds', '-S', action="store_true", default=None, dest="cumulative_s",
                  help='Generate cumulative scaffolds output (default=%default)')
parser.add_option('--agp', '-f', action="store", default=None, dest="agp",
                  help='AGP filename used to build up the assembly (default=%default)')
parser.add_option('--chart_output', '-t', action="store", default=None, dest="chart_output",
                  help='Cumulative chart output name (should end with png) (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply assembly fasta file.")
    sys.exit(1)

def __find_nvalue_point(cumulative_data, nvalue):
    count = 0
    prev = 0 
    for point in cumulative_data:
        if int(point) - prev == nvalue:
            return count, point
        else:
            count += 1
            prev =point
            

def __plot_data(a,num_rows,cumulative_data,output):
    chart = "Cumulative_Stats_Chart"

    if output:
        if not re.search('\.png',output):
            output += '.png'
        chart = output
    fig = plt.figure()
    num_columns = 1
    rows_cols = str(num_rows) + str(num_columns)
    ylabel = "Length (bp)"
    fig.subplots_adjust(hspace=.5)
    keyslist = cumulative_data.keys()
    keyslist.sort()

    for i in range(int(num_rows)):
        subplot_string = rows_cols + str(i + 1)
        plt.subplot(subplot_string)

        xlabel = "Number of " + re.sub(".*\s+","",keyslist[i])

        # Add point at 0 if it doesn't exist
        if cumulative_data[keyslist[i]][0] != 0:
            cumulative_data[keyslist[i]].insert(0,0)

        # Get N50/N90 values
        if re.search("Contig", xlabel):
            n50 = a.contigN50()
            n90 = a.contigN90()
        else :
            n50 = a.scaffoldN50()
            n90 = a.scaffoldN90()

        if n50:
            n50_x, n50_y =  __find_nvalue_point(cumulative_data[keyslist[i]], n50)
        if n90:
            n90_x, n90_y =  __find_nvalue_point(cumulative_data[keyslist[i]], n90)
        
        plt.title(keyslist[i])
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        N50_line = None
        N90_line = None
        color = "g"
        type = "solid"
        prev = 0

        for nValue in [n50, n90]: # plot lines
            n_x, n_y = __find_nvalue_point(cumulative_data[keyslist[i]], nValue)
            if nValue == prev:
                type = "dashed"

            if color=="g":
                N50_line = Line2D([n_x, 0], [n_x, n_y], color=color, ls="dashed")
            else:
                N90_line = Line2D([n_x, 0], [n_x, n_y], color=color, ls="dashed")

            plt.annotate("{:,}".format(int(nValue)), xy=(n_x, n_y), xytext=(5, -15), textcoords='offset points')

            if nValue == n_y:
                N_rect = Rectangle((0, 0), n_x, n_y, facecolor="w", edgecolor=color, ls=type)
                gca().add_patch(N_rect)
            else:
                plt.annotate("", xy=(n_x, n_y), xycoords='data',
                    xytext=(n_x, 0), textcoords='data',
                    arrowprops=dict(arrowstyle="->", color=color, ls="dashed")
                )
            prev = nValue

            color = "b" #Change color for next line


        plt.legend((N50_line, N90_line), ('N50', 'N90'), loc=(0.8, 0.05))
        plt.plot(cumulative_data[keyslist[i]],'.-', color="black")

    plt.savefig(chart,dpi=100)

def main():
    name = "No_Name"
    assembler = "None_Given"
    if options.name:
        name = options.name
    if options.assembler:
        assembler = options.assembler

    a = Assembly(args[0],name,assembler,options.minGap,options.minContig,options.minScaffold,options.agp)
    assembly = []
    assembly.append(a)
    s = AssemblyStatsUtil(assembly)

    title = "Basic Assembly Stats"
    headers = s.get_assembly_names()
    data = s.get_stats()
    st = SimpleTable(headers,data,title)
    output = None
    if options.output:
        output = options.output + ".basic_assembly_stats"
    st.print_output(output,options.html)

    cumulative_data = {}
    num_charts = 0
    if options.cumulative_c:
        key = "Cumulative Contigs"
        cumulative_data[key] = []
        cumulative_data[key] = a.get_cumulative_contigs()
        num_charts += 1
        
    if options.cumulative_s:
        key = "Cumulative Scaffolds"
        cumulative_data[key] = []
        cumulative_data[key] = a.get_cumulative_scaffolds()
        num_charts += 1
                
    if num_charts:
        __plot_data(a,num_charts,cumulative_data,options.chart_output)

    return 0
        
if __name__ == "__main__":
    sys.exit(main())



