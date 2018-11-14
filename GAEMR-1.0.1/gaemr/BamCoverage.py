#
# Class for generating bam coverage data
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.


from math import sqrt
import re
import os
import pysam
from RunCommand import RunCommand
from Feature import Feature
import PlatformConstant as pc
constant = pc.PlatformConstant()

class BamCoverage():
    """ This class represents coverage in a bam file. """

    def __init__(self, bam, prefix=None, window=None, physical_coverage = False, max_insert_size = None):
        self.output_prefix = self.__get_bam_basename(bam)
        self.bam = bam
        self.__check_for_index(bam)
        if prefix:
                self.output_prefix = prefix

        self.has_window_coverage = 0
        self.seq_coverage = {}
        self.window_size = None
        self.mean = 0
        self.median = 0
        self.mode = 0
        self.std_dev = 0
        self.physical_coverage = physical_coverage
        self.max_insert_size = constant.MAX_INSERT_SIZE
        if max_insert_size:
           self.max_insert_size = max_insert_size
        if window:
            self.seq_coverage = self.__get_window_coverage_from_bam(window,physical_coverage)
            self.has_window_coverage = 1
            self.window_size = window
        else:
            self.seq_coverage = self.__get_coverage_from_bam(physical_coverage)

    def __make_bam_index(self,bam):
        rc = RunCommand([constant.SAMTOOLS,"index",bam])
        rc.run_command()

    def __check_for_index(self,bam):
        bam_file_index1 = re.sub("$",".bai",bam)
        if not os.path.exists(bam_file_index1):
            self.__make_bam_index(bam)

    def get_bam_file(self):
        return self.bam

    def get_output_prefix(self):
        return self.output_prefix
    
    def __get_bam_basename(self,bam):
        return re.sub(".bam","",bam)

    def __build_bam_header_command(self,bam):
        return [constant.SAMTOOLS,"view","-H",bam]
    
    def __build_mpileup_command(self,bam):
        return [constant.SAMTOOLS,"mpileup",bam]

    def __get_lengths_from_headers(self,rc_output):
        lengths = {}

        for line in rc_output.stdout.readlines():
            if re.match('@SQ',line):
                tmp = line.split('\t')
                seq = re.sub("SN:","",tmp[1])
                length = re.sub("LN:","",tmp[2])
                lengths[seq] = int(length)

        return lengths

    def __initialize_windows(self,bam,window):
        window_coverage = {}
        lengths = self.__get_header_info(bam)
        
        for seq in sorted(lengths.iterkeys()):

            length = lengths[seq]

            for i in xrange(0,length,window):
                start = i
                end = i + window
                if end > length:
                    end = length
                if seq not in window_coverage:
                    window_coverage[seq] = []
                window_coverage[seq].append(Feature(start,end,0,seq))

        return window_coverage,lengths

    def __get_header_info(self,bam):
        cmd = self.__build_bam_header_command(bam)
        rc = RunCommand(cmd)
        return self.__get_lengths_from_headers(rc.run_command(0))        

    def get_median_coverage(self):
        return self.median

    def get_average_coverage(self):
        return self.mean

    def get_coverage_mode(self):
        return self.mode

    def get_coverage_stddev(self):
        return self.std_dev

    def __set_median(self,coverage_list=None,median=None):
        tmp = sorted(coverage_list)
        if tmp:
            self.median = tmp[len(tmp)/2]
        if median:
            self.median = median

    def __set_mean(self,coverage_list=None,mean=None):
        if coverage_list:
            self.mean = int(float(sum(coverage_list)) / len(coverage_list))
        if mean:
            self.mean = int(mean)

    def __set_mode(self,coverage_list=None,mode=None):
        if coverage_list:
            mode = 0
            freq = 0
            cov_dict = {}
            for i in coverage_list:
                if i not in cov_dict:
                    cov_dict[i] = 0
                cov_dict[i] += 1
                if cov_dict[i] > freq:
                    freq = cov_dict[i]
                    mode = i
            self.mode = mode
        if mode:
            self.mode = mode

    def __set_stddev(self,coverage_list=None,stddev=None):
        if coverage_list:
            mean = self.get_average_coverage()
            cov_len = len(coverage_list)
            std = 0
            for i in coverage_list:
                std += (i - mean) ** 2
            self.std_dev = "%.2f" % (sqrt(std / float(cov_len - 1)))
        if stddev:
            self.std_dev = stddev

    def __set_coverage_stats(self,coverage_list):
        self.__set_median(coverage_list)
        self.__set_mean(coverage_list)
        self.__set_mode(coverage_list)
        self.__set_stddev(coverage_list)

    def __get_sequence_coverage(self, samfile, window_coverage,lengths, window):
        total_coverage = []

        for id in sorted(lengths.iterkeys()):
            total = 0
            coverage = 0
            for pileupcolumn in samfile.pileup(id,0,lengths[id]):
                depth = int(pileupcolumn.n)
                if pileupcolumn.pos <= lengths[id]:
                    total_coverage.append(depth)
                    window_coverage[id][int(pileupcolumn.pos/window)].increment_value(depth)
                    coverage += depth
                    total += 1

        return window_coverage, total_coverage

    def __get_min_max(self,a,b):
        if a < b:
            return a,b
        return b,a

    def __get_physical_coverage(self, samfile, window_coverage,lengths, window):
        total_coverage = []
        tmp_start = 0

        #cov_dict = {}
        #size_sum = 0
        #size_total = 0
        
	for id in sorted(lengths.iterkeys()):
            seen = {}

            tmp_coverage = [0] * lengths[id]
            total_coverage += [0] * lengths[id]

            for read in samfile.fetch(id,0,lengths[id]):
                tmp = str(read).split('\t')

                if tmp[0] in seen:
                    continue
        
                if read.is_paired and read.is_proper_pair and not read.mate_is_unmapped:
                    insert_size = read.isize
                    if insert_size < self.max_insert_size:
                        start, end = self.__get_min_max(int(tmp[3]),int(tmp[7]) + int(tmp[8]))
                        if end >= lengths[id]:
                            end = lengths[id] - 1
                        #size_sum += end - start
                        #size_total += 1
                        #print "SIZE: " + str(end - start)
                        tmp_coverage[start] += 1
                        tmp_coverage[end] -= 1
                seen[tmp[0]] = 1

            cvg_sum = 0

            for x in xrange(len(tmp_coverage)):
                cvg_sum += tmp_coverage[x]
            
                total_coverage[tmp_start + x] = cvg_sum
                #if cvg_sum not in cov_dict:
                #    cov_dict[cvg_sum] = 0
                #cov_dict[cvg_sum] += 1
            
                window_coverage[id][int(x/window)].increment_value(tmp_coverage[x])
            tmp_start += lengths[id] 
        
	#print "SSUM:",
        #print size_sum,
        #print "STOT:",
        #print size_total
        
        #sum = 0
        #total = 0
        #for i in xrange(len(total_coverage)):
        #    sum += total_coverage[i]
        #    total += 1
        #print sum
        #print total

        #for j in cov_dict.keys():
        #    print j,
        #    print cov_dict[j]

        return window_coverage, total_coverage


    def __get_window_coverage_from_bam(self,window,physical_coverage):
        bam = self.get_bam_file()

        window_coverage,lengths = self.__initialize_windows(bam,window)

        samfile = pysam.Samfile(bam,"rb")

        if physical_coverage:
            seqs, total_coverage = self.__get_physical_coverage(samfile, window_coverage,lengths, window)
        else:
            seqs, total_coverage = self.__get_sequence_coverage(samfile, window_coverage,lengths, window)

        samfile.close()

        self.__set_coverage_stats(total_coverage)

        return seqs

    def __get_coverage_from_bam(self, physical_coverage):
        bam = self.get_bam_file()

        lengths = self.__get_header_info(bam)
        sorted_lengths = sorted(lengths.values(),reverse=True)
        window = sorted_lengths[0]
        return self.__get_window_coverage_from_bam(window,physical_coverage)

    def get_window_dict_from_table(self):
        windows = {}
        if self.__has_window_coverage():
            tmp = self.get_window_coverage_table()
            for i in tmp:
                if i not in windows:
                    windows[i] = {}
                for j in tmp[i]:
                    windows[i][j] = 1
        return windows

    def __has_window_coverage(self):
        return self.has_window_coverage
    
    def get_window_coverage_table(self):
        if self.__has_window_coverage:
            return self.seq_coverage
        return self.get_seq_coverage_table()

    def __get_cov_data(self):
        return self.seq_coverage

    def get_seq_coverage_table(self):
        coverage = {}
        data  = self.__get_cov_data()
        for id in sorted(data.iterkeys()):
            coverage[id] = 0
            t_cov = 0
            t_len = 0
            for i in xrange(len(data[id])):
                t_cov += data[id][i].get_value()
                t_len += data[id][i].get_length()
            coverage[id] = "%.1f" % (float(t_cov)/t_len)
        return coverage

    def get_coverage_stats(self):
        stats = []
        stats.append(["Average Cvg", self.get_average_coverage()])
        stats.append(["Cvg Std Dev", self.get_coverage_stddev()])
        stats.append(["Median Cvg", self.get_median_coverage()])
        stats.append(["Cvg Mode", self.get_coverage_mode()])
        return stats

    def __print_cov_data(self,cov_dict):
        for id in sorted(cov_dict.iterkeys()):
            print id, str(cov_dict[id])

    def print_coverage_table(self):
        self.__print_cov_data(self.get_seq_coverage_table())
