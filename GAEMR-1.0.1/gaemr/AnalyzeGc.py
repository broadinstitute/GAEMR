#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import ChartUtilities
import PlatformConstant as constant
import SimpleFastaFile as fasta_tools


class AnalyzeGc:
    """To analyze GC content related files"""

    #initialize class fields and populate them
    def __init__(self,fasta_representation,prefix=None,window=None):
        """initialize class fields"""
        self.chart_tools = ChartUtilities.ChartUtilities()
        self.fasta = fasta_representation
        if prefix:
            self.output_prefix = prefix
        self.contig_gc_table = dict()
        self.has_window_gc = 0
        self.window_gc = {}
        self.window_size = None
        if window:
            self.window_gc = self.__get_window_gc(window)
            self.has_window_gc = 1
            self.window_size = window
        self.contig_gc_table = self.__populate_contig_gc_table__()

    def __get_window_gc(self,window):
        window_gc = {}
        for contig_name,seq in self.fasta.get_sequence().iteritems():
            window_num = 0
            id = contig_name.replace("\n","")
            window_gc[id] = {}
            for i in xrange(0,len(seq),window):
                end = i + window
                if end > len(seq):
                    end = len(seq)
                window_gc[id][window_num] = 0
                window_gc[id][window_num] = self.__calculate_gc__(seq[i:end])
                window_num += 1
        return window_gc

    #populate GC table
    def __populate_contig_gc_table__(self):
        contig_gc_table = {}
        for contig_name,seq in self.fasta.get_sequence().iteritems():
            contig_gc_table[contig_name.replace("\n","")] = self.__calculate_gc__(seq)
        return contig_gc_table

    #calculate GC% of a given sequence
    def __calculate_gc__(self,seq):
        """calculate GC of a sequence"""
        to_check = seq.upper()
        gc_count = to_check.count('G') + to_check.count('C')
        return (float(gc_count)/len(seq))*100

    def __has_gc_window(self):
        return self.has_window_gc

    #public accessor to field
    def get_contig_gc_table(self):
        """return internal GC to contig table data structure"""
        return self.contig_gc_table

    def get_window_gc_table(self):
        if self.__has_gc_window():
            return self.window_gc
        return self.get_contig_gc_table()

    # private function to test if initialized with output prefix
    def __has_out_prefix(self):
        return self.output_prefix

    def get_window_dict_from_table(self):
        windows = {}
        if self.__has_gc_window():
            tmp = self.get_window_gc_table()
            for i in tmp:
                if i not in windows:
                    windows[i] = {}
                for j in tmp[i]:
                    windows[i][j] = 1
        return windows

    #public accessor to basic contig GC plot
    def get_contig_gc_plot(self):
        """Plot of graph of contig GC%"""
        labels=[]
        values=[]
        index=0
        output_file = None
        if self.__has_out_prefix():
            for contig_name,gc in self.get_contig_gc_table().iteritems():
                labels.insert(index,contig_name)
                values.insert(index,gc)
                index += 1
            output_file = self.chart_tools.gen_bar_chart_with_axis(
                values,labels,"Contig","GC%","Contig GC%",self.output_prefix + ".plot")
        return output_file

    #public accessor to basic contig GC histogram
    def get_contig_gc_histogram(self):
        """Plot of graph of contig GC%"""
        output_file = None
        if self.__has_out_prefix():
            output_file = self.chart_tools.gen_histogram(
                self.get_contig_gc_table().values(),"GC% Bin","GC% Frequency","Contig GC% Histogram",self.output_prefix + ".histo")
        return output_file








#----tests----

#fasta = fasta_tools.SimpleFastaFile("/gsap/assembly_pipeline/projects/B453/D1347/D1347_allpaths_100x50x0x_3465/run/ASSEMBLIES/test/analysis/analysis/submission.contigs.fasta")

#gc = AnalyzeGc(fasta, "analyze_gc",1000)
#print gc.get_window_gc_table()
#print gc.get_contig_gc_plot()
#print gc.get_contig_gc_histogram()
