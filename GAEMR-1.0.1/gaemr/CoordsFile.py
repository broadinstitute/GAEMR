#
# Class for manipulating nucmer coords files
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import os
import operator
from Bio import SeqIO
import PlatformConstant as pc
constant = pc.PlatformConstant()

class CoordsFile():
    """ This class represents a manipulation of nucmer coords files. """

    def __init__(self, file):
        self.ref_lengths = []
        self.query_lengths = []
        self.ref = None
        self.query = None
        self.query_gap_info = {}
        self.ref_gap_info = {}
        (self.ref, self.query, self.coords, self.ref_lookup, self.query_lookup) = self.__read_coords_file(file)

    def __set_ref_file(self, file):
        self.ref = file

    def __set_query_file(self, file):
        self.query = file

    def __set_ref_lengths(self, l_dict):
        self.ref_lengths = l_dict

    def __set_query_lengths(self, l_dict):
        self.query_lengths = l_dict
        
    def __get_base_name(self, file):
        return re.sub("\.f.*a$","",os.path.basename(file))

    def get_ref_name(self):
        return self.__get_base_name(self.ref)

    def get_query_name(self):
        return self.__get_base_name(self.query)
    
    def __get_lengths(self, file):
        l_list = []
        for seq in SeqIO.parse(open(file, 'r'), "fasta"):
            l_list.append((seq.id, len(seq)))
        return l_list

    def __format_path(self, f_tuple):
        return (os.path.abspath(f_tuple[0]),os.path.abspath(f_tuple[1]))

    def __initialize_cov_dict(self, l_list):
        c_dict = {}
        for (id, length) in l_list:
            c_dict[id] = [0] * (length + 1)
        return c_dict
        
    def __read_coords_file(self,file):
        try:                                # open up input
            in_fh  = open(file, 'r')
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            sys.exit(-1)

        ref = query = None
        coords = []
        ref_dict = {}
        query_dict = {}

        ref_cov_dict = {}
        query_cov_dict = {}

        skip_line = 1
        first = 1
        for line in in_fh.readlines():
            if re.match("^=",line):
                skip_line = 0
                continue
            if skip_line:
                if first:
                    (ref,query) = line.rstrip('\n').split()
                    self.__set_ref_file(ref)
                    self.__set_query_file(query)
                    self.__set_ref_lengths(self.__get_lengths(self.ref))
                    self.__set_query_lengths(self.__get_lengths(self.query))
                    ref_cov_dict = self.__initialize_cov_dict(self.ref_lengths)
                    query_cov_dict = self.__initialize_cov_dict(self.query_lengths)
                    first = 0 
                continue

            tmp_list = line.rstrip('\n').split()
            if int(tmp_list[3]) < int(tmp_list[4]):
                dir = 0 
            else:
                dir = 1 
                (tmp_list[3], tmp_list[4]) = (tmp_list[4], tmp_list[3])

            ref_cov_dict[tmp_list[11]][int(tmp_list[0]) - 1] += 1
            ref_cov_dict[tmp_list[11]][int(tmp_list[1])] += -1
            query_cov_dict[tmp_list[12]][int(tmp_list[3]) - 1] += 1
            query_cov_dict[tmp_list[12]][int(tmp_list[4])] += -1
            
            line_list = [int(tmp_list[0]), int(tmp_list[1]), int(tmp_list[3]), int(tmp_list[4]), int(tmp_list[6]),
                          int(tmp_list[7]), float(tmp_list[9]), tmp_list[11], tmp_list[12], dir]

            coords.append(line_list)

            if tmp_list[11] not in ref_dict:
                ref_dict[tmp_list[11]] = []
            ref_dict[tmp_list[11]].append(line_list)

            if tmp_list[12] not in query_dict:
                query_dict[tmp_list[12]] = []
            query_dict[tmp_list[12]].append(line_list)

        self.__find_gaps_from_cov(ref_cov_dict, self.ref_lengths, True)
        self.__find_gaps_from_cov(query_cov_dict, self.query_lengths, False)

        return ref, query, coords, ref_dict, query_dict

    def __find_gaps_from_cov(self, c_dict, l_list, is_ref=False):
        
        for (id, length) in l_list:
            gaps = []
            novel = []
            details = []
            
            cvg_sum = c_dict[id][0]
            start = 0
           
            in_gap = False
            if not sum:
                start = 1
                in_gap = True

            for i in xrange(1,length):
                if not cvg_sum:
                    if not in_gap:
                        start = i
                        in_gap = True
                else:
                    if in_gap:
                        size = i - start
                        gaps.append(size)
                        novel.append([id, start, i - 1, size])
                        in_gap = False
                cvg_sum += c_dict[id][i]

            if in_gap:
                size = length - start + 1
                gaps.append(size)
                novel.append([id, start, length, size])                

            total = length - sum(gaps)
            details.append([id, length, total, str("%.2f" % ((total/float(length))*100))])

            self.__update_gap_details(id, gaps, novel, details, is_ref)
            
    def get_coords(self):
        return self.coords

    def set_coords(self, coords=[]):
        self.coords = coords

    def filter_self_align(self):
            
        filtered = []
        for line in self.get_coords():
            if line[0] == line[2] and line[1] == line[3] and line[4] == line[5] and line[6] == '100.00' and line[7] == line[8]:
                continue
            filtered.append(line)
        return filtered

    def filter_by_pct_id_and_len(self, id=0, len=0):
            
        filtered = []
        for line in self.get_coords():
            if float(line[6]) > float(id) and int(line[5]) > int(len):
                filtered.append(line)
        return filtered
    
    def filter_by_pct_id(self, id=0):
        return self.filter_by_pct_id_and_len(id, '0')

    def filter_by_hit_len(self, len=0):
        return self.filter_by_pct_id_and_len('0', len)
    
    def __to_string_hits(self, t_list, delimiter=constant.TABLE_DELIMITER):
        string = ''
        for index, line in enumerate(t_list):
            if line[-1]:
                (line[2],line[3]) = (line[3], line[2])
            string += ''.join([str(x) + delimiter for x in line[:-1]])
            string = string[:-1] + '\n'
        return string[:-1]

    def __get_line_from_id(self, t_dict, id):
        if id in t_dict:
            return t_dict[id]
        return []
        
    def get_coords_by_query_id(self, query_id):
        return self.__get_line_from_id(self.query_lookup,query_id)

    def get_coords_by_ref_id(self, ref_id):
        return self.__get_line_from_id(self.ref_lookup, ref_id)

    def print_coords(self, delimiter=constant.TABLE_DELIMITER):
        return self.__to_string_hits(self.get_coords(), delimiter)

    def __sort_coords(self, item1=0, item2=0):
        coords = self.get_coords()

        if not len(coords):
            return []

        return sorted(coords,key=operator.itemgetter(item1,item2))

    def sort_by_query_start(self):
        return self.__sort_coords(8, 2)

    def sort_by_query_end(self):
        return self.__sort_coords(8, 3)

    def sort_by_ref_start(self):
        return self.__sort_coords(7, 0)
    
    def sort_by_ref_end(self):
        return self.__sort_coords(7, 1)
    
    def __get_length_by_id(self, id, is_ref=None):
        length_tuple_list = self.query_lengths
        if is_ref:
            length_tuple_list = self.ref_lengths

        for (l_id, l_length) in length_tuple_list:
            if id == l_id:
                return l_length

    def get_ref_length_by_id(self, id):
        return self.__get_length_by_id(id, True)

    def get_query_length_by_id(self, id):
        return self.__get_length_by_id(id, False)

    def get_ref_length(self, id=None):
        if id:
            return get_ref_length_by_id(id)
        return sum([x[1] for x in self.ref_lengths]) 

    def get_query_length(self, id=None):
        if id:
            return get_query_length_by_id(id)
        return sum([x[1] for x in self.query_lengths])

    def __update_gap_details(self, id, gaps, novel, details, is_ref=False):
        gap_dict = self.query_gap_info
        if is_ref:
            gap_dict = self.ref_gap_info
        if id not in gap_dict:
            gap_dict[id] = []
        gap_dict[id].append(gaps)
        gap_dict[id].append(novel)
        gap_dict[id].append(details)
        
    def get_query_gap_sizes(self):
        return self.__gap_details_helper(self.query_gap_info, self.query_lengths, 0)
    
    def get_ref_gap_sizes(self):
        return self.__gap_details_helper(self.ref_gap_info, self.ref_lengths, 0)

    def get_pct_query_uncovered(self):
        return ("%.2f" % ((sum(self.get_query_gap_sizes())/float(self.get_query_length()))*100))
        
    def get_pct_query_covered(self):
        return ("%.2f" % (100.0 - float(self.get_pct_query_uncovered())))

    def get_pct_ref_uncovered(self):
        return ("%.2f" % ((sum(self.get_ref_gap_sizes())/float(self.get_ref_length()))*100))

    def get_pct_ref_covered(self):
        return ("%.2f" % (100.0 - float(self.get_pct_ref_uncovered())))

    def __gap_details_helper(self, gap_dict, lengths_dict, list_pos):
        gaps = []

        for (id, length) in lengths_dict:
            gaps += gap_dict[id][list_pos]
        return gaps

    def get_query_coverage_details(self):
        return self.__gap_details_helper(self.query_gap_info, self.query_lengths, 2)

    def get_ref_coverage_details(self):
        return self.__gap_details_helper(self.ref_gap_info, self.ref_lengths, 2)

    def get_novel_query_details(self):
        return self.__gap_details_helper(self.query_gap_info, self.query_lengths, 1)

    def get_novel_ref_details(self):
        return self.__gap_details_helper(self.ref_gap_info, self.ref_lengths, 1)

    def __pct_id_helper(self, lengths_dict):
        pct_id = 0.0
        total = 0
        for (id, length) in lengths_dict:
            coords = self.get_coords_by_query_id(id)
            for line in coords:
                pct_id += float(re.sub("%","", str(line[6])))
                total += 1
        return pct_id/total
        
    def get_pct_identity(self):
        return ("%.2f" % (self.__pct_id_helper(self.query_lengths)))

    def get_refs_with_hits(self):
	return self.ref_lookup.keys()

    def get_querys_with_hits(self):
	return self.query_lookup.keys()
                
################## TEST ############################
#coords = "/home/unix/ssykes/dev/python/TEST.coords"
#coords = "/gsap/assembly_pipeline/projects/B472/G17679/G17679_allpaths_50f50j_4886/run/ASSEMBLIES/test/analysis/analysis/Ref_covered.coords"
#coords = "/gsap/assembly_analysis/projects/B318/pcr_free_assemblies/SEAN.coords"
#c = CoordsFile(coords)
#print c.get_ref_name()
#print c.get_query_name()
#print c.get_ref_length()
#print c.get_query_length()
#print c.get_ref_length_by_id('contig1')
#print c.get_ref_length_by_id('gi|254160123|ref|NC_012967.1|')
#print c.get_query_length_by_id('contig000009')
#print c.get_query_length_by_id('contig000100')
#c.sort_by_query_start()
#print c.print_coords()
#print c.get_query_gap_sizes()
#print c.get_ref_gap_sizes()
#print c.get_pct_query_uncovered()
#print c.get_pct_query_covered()
#print c.get_pct_ref_uncovered()
#print c.get_pct_ref_covered()
#print c.get_refs_with_hits()
