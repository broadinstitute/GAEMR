#
# Class for manipulating picard output metrics files
#

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

class PicardFile:
    """ This class represents a Picard metrics file which has certain output formats. """

    # Picard files have the characteristics as such:
    # 1. Header information lines starting with #
    #    a. contains run information
    #    b. contains date run information
    # 2. blank lines
    # 3. Metrics information if available:
    #    a. Header line
    #    b. One or more lines of metrics
    # 4. blank lines
    # 5. Histogram data
    #    a. Header line
    #    b. One or more columns; column 1 is name
    
    def __init__(self,filename):
        self.filename = filename
        self.command_line,self.date_run,self.m_headers,self.metrics,self.h_headers,self.histo = self.__load_file(filename)

    # private function to load the picard output file
    def __load_file(self,filename):
        """Tries to open picard file and starts the parser"""
        try:
            f = open(filename,'r')
            lines = f.readlines()
            f.close()

        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1

        blanks = self.__find_blank_lines(lines)
        return self.__parse_file(blanks,lines)

    # private function to find the lines which are blank in file
    def __find_blank_lines(self,lines):
        """Returns list of blank lines in file"""
        blanks = []
        count = 0
        for i in lines:
            if re.search("^$",i):
                blanks.append(count)
            count += 1
        return blanks

    # private function to parse the metrics from file
    def __parse_metrics(self,lines):
        """Returns metrics dict and header list information"""
        metrics = {}
        headers = []

        for i in xrange(len(lines)):
            if not i:
                continue
            elif i == 1:
                headers = lines[i].rstrip('\n').split('\t')
                for i in headers:
                    metrics[i] = []
            else:
                tmp = lines[i].rstrip('\n').split('\t')
                for j in range(len(tmp)):
                    metrics[headers[j]].append(tmp[j])

        return headers,metrics
    
    # private function to parse the histogram data
    def __parse_histo(self,lines):
        """Returns histogram data dict and header list information"""
        histo = {}
        headers = []

        for i in xrange(len(lines)):
            if not i:
                continue
            elif i == 1:
                headers = lines[i].rstrip('\n').split('\t')
            else:
                tmp = lines[i].rstrip('\n').split('\t')
                key = None
                for j in xrange(len(tmp)):
                    if not j:
                        histo[int(tmp[j])] = []
                        key = int(tmp[j])
                    else:
                        histo[key].append(tmp[j])
                    
        return headers,histo

    # private function to parse picard file
    def __parse_file(self,blanks,lines):
        """Gets all the instance variables for the class"""
        command_line = date_run = None
        metrics = {}
        h_headers = []
        histo = {}
        m_headers = []

        start = 0
        for i in xrange(len(blanks)):
            end = int(blanks[i])
            if not i:
                command_line,date_run = self.__parse_header(lines[start:end])
            else:
                if end - start <= 1:
                    start = end + 1
                    continue
                elif re.search("METRICS",lines[start]):
                    m_headers,metrics = self.__parse_metrics(lines[start:end])
                else:
                    h_headers,histo = self.__parse_histo(lines[start:end])
            start = end + 1
            
        return command_line,date_run,m_headers,\
               metrics,h_headers,histo

    # private function to get the header information (e.g. date run, etc.)
    def __parse_header(self,header):
        """Gets the information about how picard software was run"""
        command_line = date_run = None

        for i in range(int(len(header))):
            line = self.__clean_line(header[i])
        
            if not i or i == 2:
                continue
            elif i == 1:
                command_line = line
            else:
                date_run = line

        return command_line,date_run

    # private function to remove beginning # and strip newlines
    def __clean_line(self,line):
        """Returns cleaned comment line"""
        return str(re.sub("^#\s+","",line)).rstrip('\n')

    # private function to check if we captured data
    # not all data is in every file output
    def __has_data(self,to_check):
        """Returns whether or not the python type contains data"""
        if type(to_check).__name__ == 'dict':
            return len(to_check.keys())
        elif type(to_check).__name__ == 'list':
            return len(to_check)
        return 0

    # public function to get last bin for histogram
    def get_max_histo_bin(self):
        """Returns last entry in histogram data dict"""
        if self.__has_data(self.histo):
            return int(sorted(self.histo.iterkeys())[-1])
        return -1

    # public function to get the number of histogram bins
    def get_num_bins(self):
        """Returns the total number of bins in histogram data"""
        if self.__has_data(self.histo):
            return len(self.histo.keys())
        return -1

    # public function to get histogram
    def get_histo(self):
        """Returns histo data if there is any"""
        if self.__has_data(self.histo):
            return self.histo
        return {}

    # public function to get the histogram headers
    def get_histo_headers(self):
        """Returns histogram headers if there is any"""
        if self.__has_data(self.h_headers):
            return self.h_headers
        return []

    # public function to find out how many columns are in histo
    def num_histo_columns(self):
        """Returns column count of histogram dict if there are values."""
        if self.__has_data(self.h_headers):
            return len(self.h_headers)
        return []

    # public function to get metrics headers
    def get_metrics_headers(self):
        """Returns metrics headers if there are any"""
        if self.__has_data(self.m_headers):
            return self.m_headers
        return []

    # public function to get metrics data
    def get_metrics(self):
        """Returns metrics data if there is any."""
        if self.__has_data(self.metrics):
            return self.metrics
        return {}

    # public function to get the value of a particular metric
    def get_metric_value(self,name):
        """Returns value for metric if found."""
        for i in self.metrics.keys():
            if re.search(i,name):
                return self.metrics[i]
        return None

    # public function to get the way picard was run
    def get_command_run(self):
        """Returns the picard command"""
        return self.command_line

    # public function to get the picard date run
    def get_date_run(self):
        """Returns date picard module ran"""
        return self.date_run

    # public function to get list of histo values for a column
    def get_histo_by_column(self,column):
        histo = []
        if self.num_histo_columns() >= column:
            data = self.get_histo()
            for i in data:
                histo.append(data[i][column])
        return histo

#test = PicardFile("/home/unix/ssykes/dev/python/Fragments.scaffolds.insert_size.metrics")
#test2 = PicardFile("/seq/finishing_scratch/gaag/read_analysis_revamp/AnalyzeNGS/bin/frag.frag_reads.scaffolds.sorted.bam.qualityScoreDistribution.table.txt")
