#
# Class for manipulating agp files
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
import PlatformConstant as pc
constant = pc.PlatformConstant()



class AgpFile:
    """This class represents an agp file and its fields."""
    
    def __init__(self, agp_file=None, minGapSize=None, minConSize=None, minScaffSize=None,is_version_2=False):
        self.agp_file = agp_file
        self.minGapSize = constant.MIN_GAP_SIZE
        self.minConSize = constant.MIN_CONTIG_SIZE
        self.minScaffSize = constant.MIN_SCAFFOLD_SIZE
        self.is_version_2 = is_version_2
        if minGapSize:
            self.minGapSize = minGapSize
        if minConSize:
            self.minConSize = minConSize
        if minScaffSize:
            self.minScaffSize = minScaffSize

        self.scaffolds = {}
        self.scaffold_list = []
        self.contigs = {}
        if agp_file:
            self.__read_agp_file()

    # private function to store the agp information based on scaffold id
    def __set_scaff_record(self,scaffold_id,feature_num,line):
        """Sets up scaffolds dict to hold agp information."""
        if not scaffold_id in self.scaffolds:
            self.scaffolds[scaffold_id] = {}
            self.scaffold_list.append(scaffold_id)
        self.scaffolds[scaffold_id][feature_num] = ""
        self.scaffolds[scaffold_id][feature_num] = line
        #self.scaffold_list.append(scaffold_id)

    def __set_contig_record(self,scaffold_id,s_start,s_end,contig_id):
        self.contigs[contig_id] = {}
        self.contigs[contig_id]['scaffold_id'] = scaffold_id
        self.contigs[contig_id]['scaff_start'] = s_start
        self.contigs[contig_id]['scaff_end'] = s_end

    # private function to load agp
    def __read_agp_file(self):
        """Opens agp and gets the records."""
        agp_file = self.get_agp_file()
        try:
            f = open(agp_file,'rU')
            for lines in f.readlines():
                tmp = lines.rstrip('\n')
                data = tmp.split('\t')
                self.__set_scaff_record(data[0],int(data[3]),tmp)
                if not self.is_gap(data[0],int(data[3])):
                    self.__set_contig_record(data[0],int(data[1]),int(data[2]),data[5])
        except IOError as (errno,strerror):
                print "I/O Error({0}): {1}".format(errno,strerror)
                return -1
        return 1

    # public function to get agp filename
    def get_agp_file(self):
        """Returns agp file name."""
        return self.agp_file

    # public function to get the minimum gap size
    def get_min_gap_size(self):
        """Returns minimum gap size"""
        return int(self.minGapSize)

    # public function to get minimum scaffold size
    def get_min_scaff_size(self):
        """Returns minimum scaffold size"""
        return int(self.minScaffSize)

    # public function to get minimum contig size
    def get_min_contig_size(self):
        """Returns minimum contig size."""
        return int(self.minConSize)

    # public function to print out an agp
    def print_agp(self,output=None):
        """Prints out the agp fields."""
        if self.get_agp_file:
            outfile = sys.stdout
            if output:
                try:
                    outfile = open(output,'w')
                except IOError as (errno,strerror):
                    print "I/O Error({0}): {1}".format(errno,strerror)
                    return -1

            for i in self.scaffold_list:
                for j in sorted(self.scaffolds[i].iterkeys()):
                    outfile.write(self.scaffolds[i][j] + '\n')

    #public function to print a subset of the agp
    def get_agp_file_record_subset(self, scaffold_id, start=None, end=None):
        if self.__contains_scaffold(scaffold_id):
            if start:
                start_feature = self.__get_feature_record_by_position(scaffold_id, start)
            else:
                start_feature = 1
            if end:
                end_feature = self.__get_feature_record_by_position(scaffold_id, end)+1
            else:
                end_feature = self.__get_feature_count(scaffold_id)
            print str(start_feature) + " -- " + str(end_feature)
            for feature in range(start_feature, end_feature):
                print self.__get_feature_record(scaffold_id,feature)
            return 0
        else:
            print "Error scaffold " + scaffold_id + " not found in the AGP file"
            return -1

    # public function to get list of features based on scaffold id
    def get_agp_file_record(self,scaffold_id):
        return self.scaffolds[scaffold_id]

    # private function to count the features in a given scaffold
    def __get_feature_count(self,scaffold_id):
        return len(self.get_agp_file_record(scaffold_id))

    # public function to get the scaffolds in the agp
    def get_agp_scaffolds(self):
        """Returns list of scaffolds in agp."""
        return sorted(self.scaffolds.keys())

    # private function to get the feature record based on scaffold and
    #   feature number
    def __get_feature_record(self,scaffold_id,feature_num):
        return self.scaffolds[scaffold_id][int(feature_num)]

    # private function to get the feature record based on scaffold and
    #   position
    def __get_feature_record_by_position(self, scaffold_id, position):
        for feature_id in self.get_agp_file_record(scaffold_id):
            start = int(self.__get_feature_record(scaffold_id, feature_id).split('\t')[1])
            end = int(self.__get_feature_record(scaffold_id, feature_id).split('\t')[2])
            if start <= int(position) <= end:
                return int(feature_id)
        else:
            return 0

    # public function to get the contig id from scaffold and feature num
    def get_contig_id(self,scaffold_id,feature_num):
        """Returns contig id."""
        if self.is_gap(scaffold_id,feature_num):
            return 0
        return self.__get_feature_record(scaffold_id,feature_num).split('\t')[5]

    # private function to see if scaffold exists
    def __contains_scaffold(self,scaffold_id):
        if scaffold_id in self.get_agp_scaffolds():
            return 1
        return 0

    # public function to get scaffold length from id
    def get_scaffold_length(self, scaffold_id):
        length = 0
        if self.__contains_scaffold(scaffold_id):
            for f in self.get_agp_file_record(scaffold_id):
                length += self.get_feature_length(scaffold_id,f)
        return length

    def get_number_of_contigs_in_scaffold(self,scaffold_id):
        count = 0
        if self.__contains_scaffold(scaffold_id):
            for f in self.get_agp_file_record(scaffold_id):
                if self.is_gap(scaffold_id,f):
                    continue
                count += 1
        return count
    
    # public function to get scaffold id from contig id
    def get_scaffold_from_contig_id(self, contig_id):
        if contig_id in self.contigs:
            return self.contigs[contig_id]['scaffold_id']
        return None
        
    # public function to get the length of a feature (contig or gap)
    def get_feature_length(self,scaffold_id,feature_num):
        """Returns contig or gap length"""
        if self.is_gap(scaffold_id,feature_num):
            return int(self.__get_feature_record(scaffold_id,feature_num).split('\t')[5])
        return int(self.__get_feature_record(scaffold_id,feature_num).split('\t')[7])

    # public function to tell if record is a gap
    def is_gap(self,scaffold_id,feature_num):
        """Returns whether or not the feature is a gap"""
        line = self.__get_feature_record(scaffold_id,feature_num).split('\t')
        if line[4] in ('N', 'U'):
            return 1
        return 0

    def set_version(self,version):
        self.is_version_2 = version

    def __is_version_2(self):
        return self.is_version_2

    def get_version(self,version):
        if self.__is_version_2():
            return "2.0"
        return "1.0"

    def is_contig_reverse_complemented(self, contig_id):
        scaffold_id = self.get_scaffold_from_contig_id(contig_id)
        for i in xrange(len(self.get_agp_file_record(scaffold_id))):
            index = i + 1
            if self.get_contig_id(scaffold_id, index) == contig_id:
                if self.__get_feature_record(scaffold_id, index).split('\t')[-1] == '-':
                    return 1
                break
        return 0

    def get_coordinates_from_contig_id(self, contig_id):
        if contig_id in self.contigs:
            return self.get_scaffold_from_contig_id(contig_id),\
                   self.contigs[contig_id]['scaff_start'],self.contigs[contig_id]['scaff_end']
        return None

    def make_agp_from_assembly(self,scaffolds_list,want_version_2=True,rename=True):
        nc = 0
        ns = 0
        self.set_version(want_version_2)
        for s in scaffolds_list:
            ns += 1
            pos = 1
            element = 0
            sname = s.get_name()
            if rename:
                sname = "scaffold%05d" % (ns)
                sys.stdout.write("NEW_SCAFFOLD:\t" + sname + "\t" + "OLD_SCAFFOLD:\t" + s.get_name() + "\n")
            gaps = s.get_gaps_lengths_list()
            contigs = s.get_contig_lengths_list()
            c_ids = s.get_contig_ids()
            for l in xrange(len(contigs)):
                element += 1
                nc += 1
                cname = c_ids[l]
                if rename:
                    cname = "contig%06d" % (nc)
                    sys.stdout.write("\tNEW_CONTIG:\t" + cname + "\t" + "OLD_CONTIG:\t" + c_ids[l] + "\n")
                agpfields = [sname,pos, pos + contigs[l] - 1, element, 'W', cname, 1, contigs[l], '+']
                self.__set_scaff_record(sname,element,"\t".join([str(x) for x in agpfields]))
                self.__set_contig_record(sname,pos,pos + contigs[l] - 1,cname)
                pos += contigs[l]
                if not l == len(contigs) - 1:
                    element += 1
                    agpfields = [sname,pos, pos + gaps[l] - 1, element, 'N', gaps[l]]
                    if want_version_2:
                        agpfields += ['scaffold', 'yes', 'paired-ends']
                    else:
                        agpfields += ['fragment', 'yes', '']
                    self.__set_scaff_record(sname,element,"\t".join([str(x) for x in agpfields]))
                    pos += gaps[l]
