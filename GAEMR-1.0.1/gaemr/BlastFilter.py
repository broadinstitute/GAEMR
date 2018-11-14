#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import operator
from Bio.Blast import NCBIXML
import sys

class BlastFilter(object):
    """ This class represents a filter for a BLAST alignment."""

    def __init__(self, xml_file, want_gi_ids=True):
        self.want_gi_ids = want_gi_ids
        if xml_file:
            blast_objects = self.__get_blast_obj_from_xml(xml_file)
            #     self.id_dict = self.__build_ID_lookup(blast_objects)
            (self.blast_dict, self.nohit_dict, self.id_dict, self.query_dict) = self.__format_blast_data(blast_objects)


    # private function to parse xml to Blast object
    def __get_blast_obj_from_xml(self, xml_file):
        blast_records = NCBIXML.parse(open(xml_file))
        return blast_records

    # private function to calculate % id
    def __calculate_identity(self, identities, align_length):
        return round(((float(identities) / align_length) * 100), 2)

    # private function to make a blast m8-like line from the hit
    def __make_blast_list(self, query, hit_id, hsp):
        mismatches = hsp.align_length - hsp.identities
        pct_identity = self.__calculate_identity(hsp.identities, hsp.align_length)
        return [query, hit_id, pct_identity, hsp.align_length, mismatches, hsp.gaps, hsp.query_start,
                hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, hsp.expect, hsp.score]

    # private function to get the blast hits from xml
    def __get_hits(self, query, hit_id, hsp, hits):
        hits[query]['hits'].append(self.__make_blast_list(query, hit_id, hsp))
        return hits

    # private function to format the blast object into our data structure
    def __format_blast_data(self, blast_records):
        hits = {}
        nohit = {}
        hit_id_lookup = {}
        query_id_lookup = {}
        for blast_record in blast_records:
            query = blast_record.query.split(' ')[0]
            query_length = blast_record.query_length
            query_id = blast_record.query_id
            query_id_lookup[query] = query_id
            hits[query] = {}
            hits[query]['length'] = query_length
            hits[query]['hits'] = []
            if len(blast_record.alignments) == 0:
                nohit[query] = 1
                hits[query]['hits'].append([query, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                continue
            for alignment in blast_record.alignments:
            #   id_lookup[self.__strip_pipe(alignment.hit_id)] = alignment.hit_def
                hit_id_lookup[alignment.hit_id] = alignment.hit_def
                for hsp in alignment.hsps:
                    if self.want_gi_ids:
                        hits = self.__get_hits(query, alignment.hit_id, hsp, hits)
                    else:
                        hits = self.__get_hits(query, alignment.hit_def, hsp, hits)

        return hits, nohit, hit_id_lookup, query_id_lookup

    def update_nohit_list(self):
        obj = self.blast_dict
        list = self.nohit_dict
        for query in obj:
            if not query in list:
                if len(obj[query]['hits']) == 0:
                    list[query] = 1

    # public function to remove queries without hits from object
    def remove_queries_with_no_hits(self):
        obj = self.blast_dict
        list = self.nohit_dict
        copy = obj.copy()
        for query in copy:
            if query in list or len(copy[query]['hits']) == 0:
                del obj[query]

    # public function to remove completely redundant hits from object
    def remove_redundant_hits(self):
        obj = self.blast_dict
        list = self.nohit_dict
        for query in obj:
            if not query in list:
                #print obj[query]['hits']
                tmp_sort = sorted(obj[query]['hits'], key=operator.itemgetter(0, 7), reverse=True)
                final_sort = sorted(tmp_sort, key=operator.itemgetter(0, 6))

                start = 0
                end = 0
                obj[query]['hits'] = final_sort
                copy = obj[query]['hits'][:]
                for hit in copy:
                    if start <= hit[6] and end >= hit[7]:
                        obj[query]['hits'].remove(hit)
                        continue
                    start = hit[6]
                    end = hit[7]

    #private function to extend hits to contig edge if <= min_contig_size
    def __extend_hit(self, hit, length, min_size=200):
        if hit[6] <= min_size:
            #print "Extending " + hit[0] + " : " + str(hit[6]) + " to 0"
            hit[6] = 0

        if length - hit[7] <= min_size:
            #print "Extending " + hit[0] + " : " + str(hit[7]) + " to " + str(length)
            hit[7] = length

        return hit


    #public function to merge overlapping hits
    def merge_and_extend_hits(self, min_size=1):
        obj = self.ignore_hits_by_dict(self.nohit_dict)
        extended = 0
        #print "MIN Size = " + str(min_size)
        for query in obj:
            merged_list = []
            length = obj[query]['length']
            if len(obj[query]['hits']) > 1: #has more than one hit
                tmp_sort = sorted(obj[query]['hits'], key=operator.itemgetter(0, 7), reverse=True)
                final_sort = sorted(tmp_sort, key=operator.itemgetter(0, 6))
                obj[query]['hits'] = final_sort
                merged_hit = self.__extend_hit(obj[query]['hits'][0], length, min_size)
                #print "RAW : ",
                #print obj[query]['hits'][0]

                copy = obj[query]['hits'][1:]
                for hit in copy:
                    hit = self.__extend_hit(hit, length, min_size)
                    #print "Merged : ",
                    #print merged_hit
                    if merged_hit[7] + min_size >= hit[6]: #hits overlap
                        merged_hit[7] = hit[7]
                        #print str(self.id_dict[hit[1]]) +"\t\t"+ str(self.id_dict[merged_hit[1]])
                        if not str(self.id_dict[merged_hit[1]]).count(str(self.id_dict[hit[1]])):
                            self.id_dict[merged_hit[1]] += ", " + self.id_dict[hit[1]]
                        extended += 1
                    else:
                        merged_list.append(merged_hit)
                        merged_hit = hit
                merged_list.append(merged_hit)
                obj[query]['hits'] = merged_list
            elif len(obj[query]['hits']):
                obj[query]['hits'][0] = self.__extend_hit(obj[query]['hits'][0], length, min_size)
        print "Merged and Extended " + str(extended) + " Hits"


    #public function to filter queries from the object based on a dict
    def filter_hits_by_dict(self, in_dict, keep):
        obj = self.blast_dict
        filtered_hits = {}
        for query in obj:
            if keep:
                if query in in_dict:
                    filtered_hits[query] = obj[query]
            else:
                if not query in in_dict:
                    filtered_hits[query] = obj[query]
        return filtered_hits

    #public function to remove queries from the object based on a dict
    def ignore_hits_by_dict(self, in_dict):
        return self.filter_hits_by_dict(in_dict, False)

    #public function to select queries from the object based on a dict
    def select_hits_by_dict(self, in_dict):
        return self.filter_hits_by_dict(in_dict, True)

    # private function to determine what to return
    def __decider(self, choice=None):
        if choice == 1:
            return self.blast_dict
        else:
            return self.ignore_hits_by_dict(self.nohit_dict)

    # public function to print full blast object
    def print_object(self, unaligned=None):
        obj = self.__decider(unaligned)
        for query in obj:
            print query
            print "\tLength: " + str(obj[query]['length'])
            print "\tHits:"
            for hit in obj[query]['hits']:
                #  print "\t" + self.id_dict[self.__strip_pipe(hit[1])]
                if hit[1] in self.id_dict:
                    print "\tID:\t" + self.id_dict[hit[1]]
                print "\tAlign:" + ''.join([str(i) + "\t" for i in hit])

    #public function to print hits in the m8 format
    def get_m8(self, unaligned=None):
        out_string = ""
        hits = self.__decider(unaligned)
        for query in hits:
            for hit in hits[query]['hits']:
                if isinstance(self, VecScreenFilter) or isinstance(self, MitoFilter) or isinstance(self, gContamFilter):
                    hit[1] = self.id_dict[hit[1]]
                out_string += ''.join([str(i) + "\t" for i in hit])
                out_string = out_string[:-1] + "\n"
        return out_string[:-1]

    #public function to filter by length and pct id
    def filter_by_length_and_pct_id(self, length, pct_id):
        hits = self.__decider('False')
        filtered_hits = []
        for query in hits:
            #print "HIT:", hits[query]['hits']
            for hit in hits[query]['hits']:
                if int(hit[3]) >= int(length) and float(hit[2]) >= float(pct_id):
                    filtered_hits.append(hit)
        return filtered_hits

    #public function to filter by length
    def filter_by_length(self, length):
        return self.filter_by_length_and_pct_id(length, '0.0')

    #public function to filter by pct id
    def filter_by_pct_id(self, pct_id):
        return self.filter_by_length_and_pct_id(0, pct_id)





    #public function to filter hits by pct of contig covered
    def filter_by_pct_covered(self, pct_cvd):
        hits = self.__decider('False')
        filtered_hits = []
        for query in hits:
            bases_covered = 0
            start = 0
            end = 0
            query_length = hits[query]['length']

            tmp_sort = sorted(hits[query]['hits'], key=operator.itemgetter(0, 7), reverse=True)
            final_sort = sorted(tmp_sort, key=operator.itemgetter(0, 6))

            for hit in final_sort:
                if start <= hit[6] and end >= hit[7]:
                    continue
                else:
                    print hit[3]
                    bases_covered += hit[3]
                start = hit[6]
                end = hit[7]


            print query + " - " + str(bases_covered)
            if bases_covered / query_length <= pct_cvd:
                for hit in hits[query]['hits']:
                    filtered_hits.append(hit)
        return filtered_hits

    #public function to print list of queries with no hits
    def print_no_hit_queries(self):
        for query in self.nohit_dict:
            print query

    #public function to print hits in the SubmissionPrep format
    def get_SubmissionPrep(self):
        self.get_removal_coords(True)

    #private function returns the 0-based id for a query
    def __get_id(self, query):
        return int(self.query_dict[query].replace("Query_", "")) - 1

    #public function to print hits in the general removal format
    ## contig_name\tcontig_len\thit_start\thit_end
    def get_removal_coords(self, as_id=False):
        out_list = []
        hits = self.__decider(False)
        for query in hits:
            length = hits[query]['length']
            for hit in hits[query]['hits']:
                name = hit[0]
                if as_id:
                    name = self.__get_id(hit[0])
                out_list.append([name, length, hit[6], hit[7], 0, self.id_dict[hit[1]]])
                ### print hit
        return out_list


class rRNAFilter(BlastFilter):
    def __init__(self, xml):
        super(rRNAFilter, self).__init__(xml)


class VecScreenFilter(BlastFilter):
    def __init__(self, xml):
        super(VecScreenFilter, self).__init__(xml)

    # private class to determine the strength of the hit
    def __get_hit_strength(self, hit, query_length):
        endRange = query_length - 25

        if hit[7] <= 25 or hit[6] <= 25 or hit[7] > endRange or hit[6] > endRange:
            ### Terminal Hits
            if self.__get_score(hit) >= 24:
                return "Strong"
            elif self.__get_score(hit) >= 19:
                return "Moderate"
            elif self.__get_score(hit) >= 16:
                return "Weak"
            else:
                return "No Hit"
        else:
            ### Internal Hits
            if self.__get_score(hit) >= 30:
                return "Strong"
            elif self.__get_score(hit) >= 25:
                return "Moderate"
            elif self.__get_score(hit) >= 23:
                return "Weak"
            else:
                return "No Hit"

    # private class to calc the score of the hit
    def __get_score(self, hit):
        numid = round(hit[3] * hit[2] / 100, 0)
        return 4 * numid - 3 * hit[3] - 3 * hit[5] - 2 * hit[4]

    def filter_weak_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        #    self.remove_queries_with_no_hits()
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                strength = self.__get_hit_strength(hit, hits[query]['length'])
                if strength == "Weak" or strength == "No Hit":
                    hits[query]['hits'].remove(hit)
                    continue
                    #    self.update_nohit_list()
        self.remove_queries_with_no_hits()

    def select_Illumina_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                if not re.search("Illumina", self.id_dict[hit[1]]):
                    hits[query]['hits'].remove(hit)
                    continue
        self.remove_queries_with_no_hits()


class MitoFilter(BlastFilter):
    def __init__(self, xml, pct_mito):
        super(MitoFilter, self).__init__(xml)
        self.pct_mito = pct_mito

    def select_valid_mito_hits(self):
        #mito_dict = {}
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            query_length = hits[query]['length']
            for hit in copy:
                if not re.search("mitochondrion", self.id_dict[hit[1]])  or  float(
                    (hit[7] - hit[6]) / float(query_length)) < self.pct_mito: #or not re.search("mitochondrion",hit[1])
                    hits[query]['hits'].remove(hit)
                    continue
                covered = float((hit[7] - hit[6]) / float(query_length))
                #print self.id_dict[hit[1]] + "\t" + str(covered)
            sys.stdout.flush()
        self.remove_queries_with_no_hits()


class gContamFilter(BlastFilter):
    def __init__(self, xml):
        super(gContamFilter, self).__init__(xml)

    def ignore_mito_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                if re.search("mito", self.id_dict[hit[1]]) or re.search("mito", hit[1]):
                    hits[query]['hits'].remove(hit)
                    continue
        self.remove_queries_with_no_hits()

    def filter_hits(self):
        hits = self.ignore_hits_by_dict(self.nohit_dict)
        for query in hits:
            copy = hits[query]['hits'][:]
            for hit in copy:
                length = hit[3]
                pct_id = hit[2]
                if not (pct_id >= 98 and length >= 50) or\
                   not (pct_id >= 94 and length >= 100) or\
                   not (pct_id >= 90 and length >= 200):
                    hits[query]['hits'].remove(hit)
                    continue
        self.remove_queries_with_no_hits()


class CombinedContaminationObject(BlastFilter):
    def __init__(self, obj_list):
        super(CombinedContaminationObject, self).__init__(None)
        self.blast_dict = {}
        self.nohit_dict = {}
        self.id_dict = {}
        self.query_dict = {}

        for obj in obj_list:
            self.id_dict.update(obj.id_dict)
            self.query_dict.update(obj.query_dict)
            if not isinstance(obj, MitoFilter) and not isinstance(obj, gContamFilter) and not isinstance(obj,
                VecScreenFilter) and not isinstance(obj, rRNAFilter):
                print "ERROR: Trying to combine non BlastFilter objects"
                sys.exit(-1)
            for query in obj.blast_dict:
                if not query in self.blast_dict:
                    self.blast_dict[query] = {}
                    self.blast_dict[query]['length'] = obj.blast_dict[query]['length']
                    self.blast_dict[query]['hits'] = []
                    self.blast_dict[query]['hits'] = obj.blast_dict[query]['hits'][:]
                else:
                    if obj.blast_dict[query]['length'] == self.blast_dict[query]['length']:
                    # print "APPEND:" + query
                    # print "ORIGINAL: ",
                    # print self.blast_dict[query]['hits']
                        self.blast_dict[query]['hits'] += obj.blast_dict[query]['hits'][:]
                        # print "APPENDED: ",
                        # print self.blast_dict[query]['hits']

                    else:
                        print "ERROR:  Trying to combine Blast filter objects from different queries: " + query + "has differing lengths " + str(
                            obj[query]['length']) + " != " + str(self.blast_dict[query]['length'])
                        sys.exit(-1)
