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
from gaemr.PicardFile import PicardFile
from gaemr.SimpleTable import SimpleTable
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] insert_size_output")

parser.add_option('--output_header', '-o', action="store", default=None, type='string',dest="output",
                  help='The chart output file header (default=%default)')
parser.add_option('--metrics_header', '-m', action="store", default=None, type='string',dest="m_file",
                  help='An output header for metrics in table format (default=%default)')
parser.add_option('--direction','-d', action="store", default=None, type='string', dest="direction",
                  help='Direction of read pairs (e.g. fr,rf) (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

DELIMITER = constant.TABLE_DELIMITER

if len(args) < 1:
    parser.error("Must supply a picard output insert size file.")
    sys.exit(-1)

# private function to get box data from a single histo column
def __get_single_box(data,count):
    tmp = []
    histo = data.get_histo()

    for key in sorted(histo.iterkeys()):
        for j in range(int(histo[key][count])):
            tmp.append(int(key))
    return tmp

# private function to populate box data
def __box_data(data):
    box = []
    for i in data:
        h = data[i].get_histo_headers()
        indices = len(h) - 1
        for j in range(0,indices):
            if options.direction:
                if not re.search(options.direction + "_count",h[j + 1]):
                    continue
            box.append(__get_single_box(data[i],j))
    return box

# private function to get total number of inserts for
#     normalization of output
def __get_total_inserts(data,count):
    total = 0
    for i in data:
        total += int(data[i][count])
    return total

# private function to get single histogram data
def __get_single_histo(data,count):
    tmp = []
    histo = data.get_histo()
    total = __get_total_inserts(histo,count)

    for x in range(data.get_max_histo_bin()):
        key = x
        if key in histo:
            tmp.append((float(histo[key][count])/total)*100.0)
        else:
            tmp.append(0)
    return tmp

# private function to get the histogram data ready
def __histo_data(data,legends):
    histo = []
    h_headers = []
    for i in data:
        h = data[i].get_histo_headers()
        indices = len(h) - 1
        for j in range(0,indices):
            if options.direction:
                if not re.search(options.direction + "_count",h[j + 1]):
                    continue
            histo.append(__get_single_histo(data[i],j))
            h_headers.append(legends[i] + "_" + __get_base_header(h[j + 1]))
    return histo,h_headers

# private function to try and get some more reasonable
#     names from file
def __get_base_header(header):
    return re.sub("\_.*","",re.sub(".*\.","",header))

# private function to start getting data ready
def __format_data(data,legends):
    box = __box_data(data)
    histo,h_headers = __histo_data(data,legends)
    return box,histo,h_headers

def __set_boxplot_colors(r):
    num_boxes = len(r['boxes'])
    for key in r.keys():
        color_index = 0
        num_values = len(r[key])
        count = 0
        for x in xrange(num_values):
            if x > len(constant.PLOT_COLORS):
                color_index = 0
            r[key][x].set_color(constant.PLOT_COLORS[color_index])
            if re.match("fliers",key):
                r[key][x].set_marker("+")
            if num_values > num_boxes:
                if count == 1:
                    color_index += 1
                    count = 0
                    continue
            else:
                color_index += 1
            count += 1

# private function to plot data
def __plot_data(chart,box,histo,h_headers):
    plt.figure()
    plt.subplot(211)             # the first subplot in the first figure
    #pos = arange(int(len(h_headers)))+.5
    #plt.yticks(pos,h_headers)
    #ax = plt.gca()
    #for label in ax.get_yticklabels():
    #    label.set_fontsize(8)
    #plt.boxplot(box,0,'rs',0)
    box_dict = plt.boxplot(box,0,'rs',0)
    __set_boxplot_colors(box_dict)
    ax = plt.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    plt.subplot(212)
    plt.xlabel('Insert Size (bp)')
    plt.ylabel('Pct of Total Pairs')
    color_index = 0
    for i in range(len(histo)):
        if color_index > len(constant.PLOT_COLORS):
            color_index = 0
        plt.plot(histo[i], label=h_headers[i],color=constant.PLOT_COLORS[color_index])
        color_index += 1
    plt.legend(h_headers,bbox_to_anchor=(0, 1.05, 1, .25), loc=3, borderaxespad=0.,prop={"size":10},\
                                        fancybox=True, ncol=2)
    #plt.legend(h_headers,loc=2,prop={"size":8}).draw_frame(False)
    plt.savefig(chart,dpi=100)

# strip out insert size stuff from file name
def __get_base_name(file):
    return re.sub(".*\/","",str((re.sub("\..*\.insert.*$","",file))))

def __get_header_line(data,labels):
    headers = ['ID']
    count = 0

    for j in data:
        metrics = data[j].get_metrics()
        tmp = metrics.keys()
        for i in xrange(int(len(metrics[tmp[0]]))):
            headers.append(labels[count])
        count += 1

    return headers

# public function to print metrics from file
def print_metrics(data,labels,m_output):
    headers = __get_header_line(data,labels)
    title = 'Insert Size Metrics'
    table_data = []
    metrics = data[0].get_metrics_headers()

    ordered_metrics = ['PAIR_ORIENTATION','READ_PAIRS','MEAN_INSERT_SIZE','STANDARD_DEVIATION',
                       'MEDIAN_INSERT_SIZE','MEDIAN_ABSOLUTE_DEVIATION','MIN_INSERT_SIZE',
                       'MAX_INSERT_SIZE','WIDTH_OF_50_PERCENT', 'WIDTH_OF_90_PERCENT',
                       'SAMPLE','LIBRARY','READ_GROUP']
    if len(metrics):
        for i in ordered_metrics:
            tmp = [i]
            for j in data:
                for k in data[j].get_metric_value(i):
                    tmp.append(k)
            table_data.append(tmp)

    st = SimpleTable(headers,table_data,title)
    st.print_output(m_output,options.html)

def __debug_print_list(type, to_print_list):
    for i in xrange(len(to_print_list)):
        for j in to_print_list[i]:
            print str(type),str(j)

def main():
    file_count = 0
    data = {}
    chart = "insert_size"
    legends = []
    
    if options.output:
        chart = options.output + ".insert_size.png"

    # get plot names data ready for each file
    for file in args:
        try:
            data[file_count] = {}
            p = PicardFile(file)
            data[file_count] = p
            legends.append(__get_base_name(file))
            file_count += 1
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1

    # get metrics file ready for output
    output = None
    if options.m_file:
        output = options.m_file + ".insert_size"
    print_metrics(data,legends,output)

    # get plot data
    box,histo,h_headers = __format_data(data,legends)

    #__debug_print_list("BOX",box)
    #__debug_print_list("HISTO",histo)
    #__debug_print_list("H_HEADERS",h_headers)
    __plot_data(chart,box,histo,h_headers)

    return 0

if __name__ == "__main__":
    sys.exit(main())

