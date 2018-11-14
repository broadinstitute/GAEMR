#
# Class for storing flagstat output
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
from RunCommand import RunCommand
import PlatformConstant as pc
constant = pc.PlatformConstant()
from SimpleTable import SimpleTable

class BamFlagstat():
    """This class represents data in samtools flagstat output"""

    def __init__(self,bam):
        self.bam = bam
        self.command = RunCommand(self.__build_command(bam))
        self.flagstat_output = self.command.run_command(0).stdout.readlines()
        self.stats = self.__populate_stats()

    # private function to make flagstat command
    def __build_command(self,bam):
        return [constant.SAMTOOLS,"flagstat",bam]

    # private function to parse output
    def __populate_stats(self):
        data = {}
        qc_fail = 0
        for i in self.flagstat_output:           
            tmp = i.rstrip('\n').split(' ')
            if re.match('\+',tmp[1]):
                qc_fail = tmp[2]
            qc_pass = tmp[0]
            key = ""
            for j in range(3,int(len(tmp))):
                if re.match('\(',tmp[j]):
                    if not re.search('mapQ',tmp[j]):
                        break
                    key += tmp[j] + "_"
                else:
                    key += tmp[j] + "_"
            key = key[:-1]
            data[key] = ()
            data[key] = (int(qc_pass), int(qc_fail))
        return data

    # private function to get the value from a particular stat
    def __get_value_from_key(self,key,key_tuple):
        sum = 0
        for i in key_tuple:
            sum += self.stats[key][i]
        return sum

    # private function to tell whether or not to look at
    # qc_pass, qc_fail or both (default)
    def __get_data_to_check(self,qc_pass=None):
        if qc_pass == None:
            return (0,1)
        if qc_pass == 1:
            return (0,)
        return (1,)

    # private function to help get the right stat
    def __check_helper(self,key, qc_pass=None):
        return self.__get_value_from_key(key,self.__get_data_to_check(qc_pass))

    # public functions follow to get number value from each key
    def num_total_reads(self,qc_pass=None):
        key = 'in_total'
        return self.__check_helper(key,qc_pass)

    def num_paired_in_sequencing(self,qc_pass=None):
        key = 'paired_in_sequencing'
        return self.__check_helper(key,qc_pass)

    def num_properly_paired(self,qc_pass=None):
        key = 'properly_paired'
        return self.__check_helper(key,qc_pass)

    def num_duplicates(self,qc_pass=None):
        key = 'duplicates'
        return self.__check_helper(key,qc_pass)        

    def num_read1_seqs(self,qc_pass=None):
        key = 'read1'
        return self.__check_helper(key,qc_pass)

    def num_read2_seqs(self,qc_pass=None):
        key = 'read2'
        return self.__check_helper(key,qc_pass)

    def num_singletons(self,qc_pass=None):
        key = 'singletons'
        return self.__check_helper(key,qc_pass)    

    def num_mapped(self,qc_pass=None):
        key = 'mapped'
        return self.__check_helper(key,qc_pass)

    def num_mapped_with_mate(self,qc_pass=None):
        key = 'with_itself_and_mate_mapped'
        return self.__check_helper(key,qc_pass)

    def num_mates_mapped_to_diff_chr(self,qc_pass=None):
        key = 'with_mate_mapped_to_a_different_chr'
        return self.__check_helper(key,qc_pass)

    def num_mates_mapped_to_diff_chr_qual(self,qc_pass=None):
        key = 'with_mate_mapped_to_a_different_chr_(mapQ>=5)'
        return self.__check_helper(key,qc_pass)

    # public functions follow to get common %
    def pct_mapped(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            if self.num_total_reads(qc_pass) != 0:
                return "%.2f" % ((float(self.num_mapped(qc_pass))/self.num_total_reads(qc_pass)) * 100)
        return 0
    
    def pct_properly_paired(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            if self.num_mapped(qc_pass) != 0:
                return "%.2f" % ((float(self.num_properly_paired(qc_pass))/self.num_mapped(qc_pass)) * 100)
        return 0

    def pct_aligned_with_pair(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            if self.num_mapped(qc_pass) != 0:
                return "%.2f" % ((float(self.num_mapped_with_mate(qc_pass))/self.num_mapped(qc_pass)) * 100)
        return 0

    def pct_singletons(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            if self.num_mapped(qc_pass) != 0:
                return "%.2f" % ((float(self.num_singletons(qc_pass))/self.num_mapped(qc_pass)) * 100)
        return 0

    def pct_paired_in_sequencing(self,qc_pass=None):
        if self.num_total_reads(qc_pass) != 0:
            return "%.2f" % ((float(self.num_paired_in_sequencing(qc_pass))/self.num_total_reads(qc_pass)) * 100)
        return 0

    def pct_duplicates(self,qc_pass=None):
        if self.num_total_reads(qc_pass) != 0:
            return "%.2f" % ((float(self.num_duplicates(qc_pass))/self.num_total_reads(qc_pass)) * 100)
        return 0
    
    def pct_chimeras(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            if self.num_mapped(qc_pass):
                return "%.2f" % ((float(self.num_mates_mapped_to_diff_chr(qc_pass))/self.num_mapped(qc_pass)) * 100)
        return 0

    def pct_chimeras_qual(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            if self.num_mapped(qc_pass) != 0:
                return "%.2f" % ((float(self.num_mates_mapped_to_diff_chr_qual(qc_pass))/self.num_mapped(qc_pass)) * 100)
        return 0

    # public function to get the stats keys; may be useful for user to peruse.
    def get_keys(self):
        return self.stats.keys()

    # public function to see if bam is mapped or not (prevents division by 0)
    def is_aligned(self,qc_pass=None):
        key = 'mapped'
        return self.__check_helper(key,qc_pass) != 0

    def is_paired(self):
        if float(self.pct_paired_in_sequencing()) > 0:
            return 1
        else:
            return 0

    def get_bam_name(self):
        return self.bam

    def __parenthesize(self, value):
        return ' (' + str(value) + '%)'

    def get_aligned_stats_table(self,qc_pass=None):
        if self.is_aligned(qc_pass):
            aligned_data = self.get_unaligned_stats_table(qc_pass)
            aligned_data.append(["Mapped" , str(self.num_mapped(qc_pass)) +
                                                self.__parenthesize(self.pct_mapped(qc_pass))])
            aligned_data.append(["Singletons" , str(self.num_singletons(qc_pass)) +
                                                    self.__parenthesize(self.pct_singletons(qc_pass))])
            aligned_data.append(["Mapped w/ Mate" , str(self.num_mapped_with_mate(qc_pass)) +
                                                        self.__parenthesize(self.pct_aligned_with_pair(qc_pass))])
            aligned_data.append(["Properly Paired" , str(self.num_properly_paired(qc_pass)) +
                                                         self.__parenthesize(self.pct_properly_paired(qc_pass))])
            aligned_data.append(["Cross-chromosome" , str(self.num_mates_mapped_to_diff_chr(qc_pass)) +
                                                          self.__parenthesize(self.pct_chimeras(qc_pass))])
            aligned_data.append(["Cross-chromosome (MQ >= 5)" , str(self.num_mates_mapped_to_diff_chr_qual(qc_pass)) +
                                                                    self.__parenthesize(self.pct_chimeras_qual(qc_pass))])
            return aligned_data
        return self.get_unaligned_stats_table(qc_pass)

    def get_unaligned_stats_table(self,qc_pass=None):
        unaligned_data = []
        unaligned_data.append(["Total Reads", self.num_total_reads(qc_pass)])
        unaligned_data.append(["Paired Reads", str(self.num_paired_in_sequencing(qc_pass)) +
                                                   self.__parenthesize(self.pct_paired_in_sequencing(qc_pass))])
        unaligned_data.append(["Duplicates",str(self.num_duplicates(qc_pass)) +
                                                self.__parenthesize(self.pct_duplicates(qc_pass))])
        unaligned_data.append(["Total Read 1",self.num_read1_seqs(qc_pass)])
        unaligned_data.append(["Total Read 2",self.num_read2_seqs(qc_pass)])
        return unaligned_data

