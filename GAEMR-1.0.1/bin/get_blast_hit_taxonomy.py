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
import operator
from Bio import SeqIO
from gaemr.RunCommand import RunCommand
import gaemr.PlatformConstant as pc
from gaemr.SimpleTable import SimpleTable
from gaemr.Taxonomy import Taxonomy

constant = pc.PlatformConstant()

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] <blast tab delimited> <query fasta file>")

parser.add_option('--nodes', '-n', action="store", default=constant.BLAST_NODES, type='string', dest="nodes_db",
                  help='NCBI nodes db file (default=%default)')
parser.add_option('--names', '-a', action="store", default=constant.BLAST_NAMES, type='string', dest="names_db",
                  help='NCBI names db file (default=%default)')
parser.add_option('--db', '-d', action="store", default=constant.BLAST_NT, type='string', dest="db",
                  help='NCBI nt db local location (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string',dest="output",
                   help='Output file header for contig taxonomic information (default=%default)')
parser.add_option('--heatmap_data','-m',action="store", default=None, type='string', dest="hm_data",
                  help='Output file name for heatmap-like output with contigs colored according to organism (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 2: 
    parser.error("Must supply a blast space delimited file and query fasta file.")
    sys.exit(-1)

# private function to get lengths of queries
def __get_lengths(file):
    """Returns length of blast query sequences."""
    lengths = {}
    try:
        for s in SeqIO.parse(open(file, 'r'), "fasta"):
            lengths[s.id] = len(s.seq)
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1
    return lengths

# private function to get blast information
def __get_blast_data(file):
    """Returns blast dictionary of lists based on query id."""
    data = {}
    try:
        f=open(file)
        for lines in f.readlines():
            blast_line = lines.rstrip('\n').split('\t')
            blast_line[6] = int(blast_line[6])
            assert len(blast_line) == 12, 'Must supply tab-delimited BLAST output.'
            if blast_line[0] in data:
                data[blast_line[0]].append(blast_line)
            else:
                data[blast_line[0]] = []
                data[blast_line[0]].append(blast_line)                
        f.close()        
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        return -1
    return data

# private function to get the gi number from the blast output
def __get_gis_from_blast_line(data):
    """Returns dict of gi numbers seen along with their original string."""
    gi_nums = {}
    for i in data:
        for j in data[i]:
            if re.search('\|',j[1]):
                ids = j[1].split('|')
                gi_nums[ids[1]] = j[1]
    return gi_nums

# private function to get the overlapping hits and
#   concatenate them down
def __parse_query_record(data):
    tmp = {}

    prev_end = 0
    prev_start = 0
    prev_gi = None

    
    data.sort(key=operator.itemgetter(1,int(6)))

    # go through the query's hits
    for j in data:

        if not re.search('\|',j[1]):
            continue
        
        # get some vitals about the hit
        gi = (j[1].split('|'))[1]
        start = int(j[6])
        end = int(j[7])
        length = (end - start) + 1

        # have we seen this gi before
        if prev_gi:

            # is it the same one?
            if gi == prev_gi:

                # check redundancy, just in case
                if start > prev_start and end < prev_end:
                    continue

                # add the new sequence hit to what we already have
                if start < prev_end:
                    tmp[gi] += (end - prev_end) + 1
                    prev_end = end
                else: # we just add in this length if no overlap
                    tmp[gi] += length
                    prev_start = start
                    prev_end = end

            # first time for gi, store it
            else:
                tmp[gi] = length
                prev_start = start
                prev_end = end
            prev_gi = gi
        # we don't have a previous gi
        else:
            tmp[gi] = length
            prev_gi = gi
            prev_start = start
            prev_end = end

    return tmp

# private function to get longest hit and gi
def __get_longest_covered_hit(data):
    """Returns dict of blast hits and hit lengths"""
    blast_hit = {}
    hit_len = {}

    # for each query
    for i in data:
        # get the blast hits
        tmp = __parse_query_record(data[i])
        
        hit_len[i] = {}

        # go through hits and make the hit length values, reverse sort sizes
        if tmp:
            for k in tmp:
                hit_len[i][k] = tmp[k]
            blast_hit[i] = sorted(tmp.iteritems(), key=operator.itemgetter(1), reverse=True)[0]
    return blast_hit,hit_len

# private function to get the output string for hit covered, etc.
def __get_coverage_string(lengths, blast_hit, tax_obj, gi_tax_dict, gi_nums):
    """Returns blast coverage line."""
    covered = []
    for i in sorted(lengths.iteritems(), key=operator.itemgetter(0)):
        q_id = i[0]
        q_length = i[1]

        tmp = [q_id,q_length]
        if q_id in blast_hit:
            covered_len = blast_hit[q_id][1]
            if covered_len > q_length:
                covered_len = q_length
            covered_pct = "%.2f" % ((float(covered_len)/q_length)*100)
            tmp += [covered_len, covered_pct]
        else:
            tmp += ["0","0.0"]

        if tax_obj.have_nodes() and tax_obj.have_names() and q_id in blast_hit:
            tmp.append(tax_obj.get_taxonomy_string(gi_tax_dict[blast_hit[q_id][0]]))
        elif q_id in gi_nums:
            tmp.append(gi_nums[blast_hit[q_id][0]])
        else:
            tmp.append("domain=NO_HIT")

        covered.append(tmp)

    return covered

# public function to check for presence of dbs
def check_dbs(nodes, names, db):
    """Checks whether or not we have some dbs to use for our analysis."""
    if not nodes and not names or nodes and not names or names and not nodes:
        if not constant.BLAST_NODES or not constant.BLAST_NAMES:
            print "If using local BLAST taxonomy dump files, must have both nodes and names files as input."
            sys.exit(-1)
        nodes = constant.BLAST_NODES
        names = constant.BLAST_NAMES
        
    if not db:
        if not constant.BLAST_NT:
            print "Must supply a path to local BLAST nt database."
            sys.exit(-1)
        db = constant.BLAST_NT

    if options.hm_data and not names or not nodes:
        print "Need nodes and names file for heatmap data."
        sys.exit(-1)

    return nodes, names, db

# private function to get a heatmap-like output
def __get_heatmap_string(contig, genus_dict, lengths):
    """Returns heatmap-like string for painting a contig to taxonomic value."""
    tmp_string = "#KEY:  "
    for k,v in sorted(genus_dict.items(), key=operator.itemgetter(1)):
        tmp_string += "%s=%d;" % (k,v)
    key_string = tmp_string[:-1]
    key_string += "\n"
    
    h_string = ""
    for i in sorted(contig.iterkeys()):
        h_string += "%s\t%s\t" % (i, lengths[i])
        for j in xrange(int(len(contig[i]))):
            h_string += "%.2f\t" % float(contig[i][j])
        h_string.rstrip('\t')
        h_string += "\n"
    return key_string + h_string

# private function to get the heatmap data
def __get_heatmap_data(gi_tax_dict, tax_obj, lengths, data, hit_lengths):
    contig = {}
    genus_dict = {}
    genus_dict["No Hit"] = 0
    genus_count = 1

    for i in data:

        contig[i] = [0] * 100
        tmp = {}
        for j in data[i]:
            if re.search('\|',j[1]):
                gi = (j[1].split('|'))[1]
                if gi not in tmp:
                    tmp[gi] = []
                tmp[gi].append([int(j[6]),int(j[7]),int(re.sub("\.\d+","",j[2]))])
        hits = sorted(hit_lengths[i].iteritems(), key=operator.itemgetter(1), reverse=True)
        gis = []

        for j in hits:
            gis.append(j[0])

        for k in gis:
            for t_list in tmp[k]:
                if k not in gi_tax_dict:
                    continue
                genus = tax_obj.get_genus(gi_tax_dict[k])
                if genus == 'root':
                    genus = tax_obj.get_species(gi_tax_dict[k])
                if genus not in genus_dict:
                    genus_dict[genus] = genus_count
                    genus_count += 1

                start = int((float(t_list[0])/lengths[i])*100)
                end = int((float(t_list[1])/lengths[i])*100)
                id = "00"
                if not re.match("100",str(t_list[2])):
                    id = t_list[2]

                for num in xrange(start,end):
                    if contig[i][num] == 0:
                        contig[i][num] = str(genus_dict[genus]) + "." + str(id)

    return __get_heatmap_string(contig,genus_dict,lengths)
    
def main():

    nodes_db, names_db, blast_db = check_dbs(options.nodes_db,options.names_db,options.db)
    
    data = __get_blast_data(args[0])
    lengths = __get_lengths(args[1])

    tax_obj = Taxonomy(nodes_db=nodes_db, names_db=names_db, blast_db=blast_db)
    gi_nums = __get_gis_from_blast_line(data)
    
    blast_hit,hit_lengths = __get_longest_covered_hit(data)

    gi_tax_dict = tax_obj.get_gi_tax_lookup(gi_nums)

    title = "Taxonomic Classification of BLAST Hits"
    headers = ["QueryId","QueryLen","QueryHitLen","PctCovered","TaxonomicString"]
    t_data = __get_coverage_string(lengths, blast_hit, tax_obj, gi_tax_dict, gi_nums)

    st = SimpleTable(headers,t_data,title)

    output = None
    if options.output:
        output = options.output + ".blast_hit_taxonomy"
    st.print_output(output,options.html)

    if options.hm_data:
        try:
            output = open(options.hm_data,'w')
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1            
        output.write(__get_heatmap_data(gi_tax_dict, tax_obj, lengths, data, hit_lengths))
        output.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())

