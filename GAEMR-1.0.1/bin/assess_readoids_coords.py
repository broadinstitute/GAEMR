#! /usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
from gaemr.SimpleTable import SimpleTable

parser = OptionParser(usage="usage: %prog [options]")

parser.add_option('--coords', '-c', action="store", default=None, type='string',dest="coords_file",
                  help='Nucmer coords file (default=%default)')
parser.add_option('--pairs', '-p', action="store", default=None, type='string',dest="pairs",
                  help='Paired readoid fasta file (default=%default)')
parser.add_option('--fasta', '-f', action="store", default=None, type='string',dest="fasta",
                  help='Scaffolds fasta file (default=%default)')
parser.add_option('--ref', '-r', action="store", default=None, type='string',dest="ref",
                  help='Reference fasta file (default=%default)')
parser.add_option('--length', '-l', action="store", default="100", type='int',dest="rd_ln",
                  help='Readoid length (default=%default)')
parser.add_option('--distance', '-d', action="store", default="100000", type='int',dest="distance",
                  help='Distance between readoids (default=%default)')
parser.add_option('--error', '-e', action="store", default=".05", type='float',dest="error",
                  help='Allowable deviation of distance between readoids (default=%default)')
parser.add_option('--verbose', '-v', action="store", default=None, dest="verbose",
                  help='Verbose output prefix (default=%default)')
parser.add_option('--output', '-o', action="store", default="ScaffoldAssessment", type='string', dest="output",
                  help='Output prefix (default=%default)')

(options, args) = parser.parse_args()

# check inputs or die

if not options.coords_file or not options.fasta or not options.ref or not options.pairs:
    parser.error("Must supply nucmer coords file using -c, scaffolds fasta file using -f, reference fasta file using -r," +
                " and paired read-oid fasta using -p")
    sys.exit(-1)

def readScaffolds(input):
    lengths = {}
    in_fh  = open(input, 'r')
    for record in SeqIO.parse(in_fh,'fasta'):
        seq = record.seq
        id  = record.id
        ln  = len(seq)
        lengths[id] = ln
    return lengths

def readRef(input):
    ref_lengths = {}
    in_fh  = open(input, 'r')
    for record in SeqIO.parse(in_fh,'fasta'):
        seq = record.seq
        id  = record.id
        ln  = len(seq)
        ref_lengths[id] = ln
    return ref_lengths

def readPairs(input):
    pairs = {}
    in_fh  = open(input, 'r')
    for record in SeqIO.parse(in_fh,'fasta'):
        id  = record.id
        pair_id  =  id
        pair_id  =  int(re.sub("^p","",re.sub("_.*","",pair_id)))
        if pair_id in pairs:
            pairs[pair_id]['read2'] = id
        else:
            pairs[pair_id] = {}
            pairs[pair_id]['read1'] = id
    return pairs

def readCoords(input):
    coords = {}

    try: 				# open up input
        in_fh  = open(input, 'r')
    except IOError as (errno,strerror):
        print "I/O Error({0}): {1}".format(errno,strerror)
        sys.exit(-1)

    skip_line = 1
    for line in in_fh.readlines():
        if re.match("^=",line):
            skip_line = 0
            continue
        if skip_line:
            continue

        tmp_list =    line.rstrip('\n').split()
        read_id  =    tmp_list[12]
        ref      =    tmp_list[11]
        ref_st   =    tmp_list[0]
        ref_end  =    tmp_list[1]
        read_st  =    tmp_list[3]
        read_end =    tmp_list[4]
        pair_id  =    read_id
        pair_id =     int(re.sub("p","",re.sub("_.*","",pair_id)))

        if read_id in coords:
            coords[read_id]['multimap'] = 1
            continue
        else:
            coords[read_id] = {}
            coords[read_id]['multimap'] = 0
            if read_st < read_end:
                coords[read_id]['dir']     = "plus"
                coords[read_id]['ref']     =  ref
                coords[read_id]['ref_pos'] =  int(ref_st)
            else:
                coords[read_id]['dir']     = "minus"
                coords[read_id]['ref']     =  ref
                coords[read_id]['ref_pos'] =  int(ref_end)
    return coords


def get_filehandle(output):
    output = re.sub("$",".pair_class_details",output)
    out_fh = output
    if output:
        try:
            out_fh = open(output,'w')
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            sys.exit(-1)
    return out_fh

def assessPairs(pair_dict,coords_dict,distance,error,scaf_lengths,ref_lengths,output,verbose):
    max_dist = distance + (distance*error)
    min_dist = distance - (distance*error)
    pair_count, unaligned_count, multimap_count, cross_chrom_count, same_orientation   = 0, 0, 0, 0, 0
    valid_internal, invalid_internal, valid_circular, invalid_orientation, other_count = 0, 0, 0, 0, 0
    if verbose:
        out_fh = get_filehandle(options.verbose)
    for pair in sorted(pair_dict.iterkeys()):
        pair_count += 1
        r1, r2  = pair_dict[pair]['read1'], pair_dict[pair]['read2']
        tmp     = r1.split("_")						# scaffold_id contained in read name:  p451_scaffold00001_4610000
        scaf    = tmp[1]
        scaf_ln = scaf_lengths[scaf]
        pair_class = " - "
        dist       = " - "

        if r1 in coords_dict:
            r1_mm  = coords_dict[r1]['multimap']
            r1_dir = coords_dict[r1]['dir']
            r1_pos = coords_dict[r1]['ref_pos']
            r1_id  = coords_dict[r1]['ref']
            r1_ref_ln  = ref_lengths[r1_id]
        else:
            r1_mm = r1_dir = r1_pos = r1_id = r1_ref_ln = '-'

        if r2 in coords_dict:
            r2_mm  = coords_dict[r2]['multimap']
            r2_dir = coords_dict[r2]['dir']
            r2_pos = coords_dict[r2]['ref_pos']
            r2_id  = coords_dict[r2]['ref']
            r2_ref_ln  = ref_lengths[r2_id]
        else:
            r2_mm = r2_dir = r2_pos = r2_id = r2_ref_ln = '-'

        #############################  Case 1:  One or both reads in pair do not align    #####################
        if r1 not in coords_dict or r2 not in coords_dict:
            unaligned_count += 1
            pair_class = "Unaligned"
        ########################  Case 2:  One or both reads in pair have multiple mapping    #################
        elif r1_mm == 1 or r2_mm == 1:
            pair_class = "MultipleMap"
            multimap_count += 1
        #####################  Case 3:  Both reads align, but to 2 different reference entries   ##############
        elif r1_id != r2_id:
            pair_class = "Invalid_CrossChromosome"
            cross_chrom_count += 1
        #################  Case 4:  Both reads align to same ref entry, but same orientation   ################
        elif r1_dir == r2_dir:
            pair_class = "Invalid_Orientation"
            invalid_orientation += 1
        #################  Case 5a:  Both reads align in valid direction (r1+,r2-) and non-circular #############
        elif r1_dir =='plus' and r1_pos < r2_pos:
            dist = (r2_pos - r1_pos) + 1
            if dist > min_dist and dist < max_dist:
                pair_class = "Valid"
                valid_internal += 1
            elif (dist > max_dist):
                pair_class = "Invalid_TooLarge"
                valid_internal += 1
            elif (dist <  min_dist):
                pair_class = "Invalid_TooSmall"
                valid_internal += 1
        #################  Case 5b:  Both reads align in valid direction (r2+,r1-) and non-circular #############
        elif r2_dir =='plus' and r2_pos < r1_pos:
            dist = (r1_pos - r2_pos) + 1
            if dist > min_dist and dist < max_dist:
                pair_class = "Valid"
                valid_internal += 1
            elif (dist > max_dist):
                pair_class = "Invalid_TooLarge"
                valid_internal += 1
            elif (dist <  min_dist):
                pair_class = "Invalid_TooSmall"
                valid_internal += 1
        #########  Case 6a:  Both reads align in invalid direction (r1+,r2-): either invalid or Circular  #########
        elif r1_dir =='plus' and r1_pos > r2_pos:
            r1_from_end = r1_ref_ln - r1_pos
            r2_from_st  = r2_pos
            dist        = r2_from_st + r1_from_end
            if (dist > min_dist and dist < max_dist):
                pair_class = "Valid_Circular";
                valid_circular += 1
            else:
                pair_class = "Invalid_Orientation"
                invalid_orientation += 1
        elif r2_dir =='plus' and r2_pos > r1_pos:
            r2_from_end = r2_ref_ln - r2_pos
            r1_from_st  = r1_pos
            dist        = r1_from_st + r2_from_end
            if (dist > min_dist and dist < max_dist):
                pair_class = "Valid_Circular";
                valid_circular += 1
            else:
                pair_class = "Invalid_Orientation"
                invalid_orientation += 1
        else:
            pair_class  = "Other"
            other_count += 1

        list = (pair, r1, r1_id, r1_pos, r1_ref_ln, r1_dir, r2, r2_id, r1_pos, r2_dir, r2_ref_ln, dist, pair_class)
        if verbose:
            for item in list:
                out_fh.write("%s\t" % item)
            out_fh.write("\n")

        #print pair, r1, r1_id, r1_pos, r1_ref_ln, r1_dir, r2, r2_id, r1_pos, r2_dir, r2_ref_ln, dist, pair_class

    pair_classes = [unaligned_count, multimap_count, cross_chrom_count, valid_internal, invalid_internal, valid_circular, invalid_orientation, other_count]
    return pair_classes,pair_count

def reportStats(list_of_classes,pair_count,output):
    title   = "Scaffold Accuracy Stats"
    headers = ['Stat', 'Total', 'Pct']
    data    = [];
    unaligned, multimap, cross_chrom, valid_internal, invalid_ln, valid_circular, invalid_orientation, other = list_of_classes

    total_valid    = valid_internal + valid_circular
    total_invalid  = cross_chrom    + invalid_ln + invalid_orientation
    pct_unaligned  = "%.3f" % ((float(unaligned)/pair_count)*100)
    pct_mm         = "%.3f" % ((float(multimap)/pair_count)*100)
    pct_cc         = "%.3f" % ((float(cross_chrom)/pair_count)*100)
    pct_invalid_ln = "%.3f" % ((float(invalid_ln)/pair_count)*100)
    pct_invalid_o  = "%.3f" % ((float(invalid_orientation)/pair_count)*100)
    pct_valid      = "%.3f" % ((float(valid_internal)/pair_count)*100)
    pct_valid_circ = "%.3f" % ((float(valid_circular)/pair_count)*100)
    pct_other      = "%.3f" % ((float(other)/pair_count)*100)
    total          = unaligned + multimap + cross_chrom + valid_internal + invalid_ln + valid_circular + invalid_orientation + other

    #print pair_count, unaligned, pct_unaligned, multimap, pct_mm, cross_chrom, pct_cc, invalid_orientation, pct_invalid_o, invalid_ln, pct_invalid_ln, valid_internal, pct_valid, valid_circular, pct_valid_circ, other, pct_other

    data.append(["Total Input Pairs", pair_count, "100.0"])
    data.append(["Multiply Mapped", multimap, pct_mm])
    data.append(["Unaligned", unaligned, pct_unaligned])
    data.append(["Cross Chromosome", cross_chrom, pct_cc])
    data.append(["Invalid Orientation", invalid_orientation, pct_invalid_o])
    data.append(["Invalid Length", invalid_ln, pct_invalid_ln])
    data.append(["Valid", valid_internal, pct_valid])
    data.append(["Valid Circular", valid_circular, pct_valid_circ])
    data.append(["Other", other, pct_other])

    st =  SimpleTable(headers,data,title)
    st.print_output(output,1)

def main():

    scaffold_length_dict  = readScaffolds(options.fasta)

    ref_length_dict  = readRef(options.ref)

    pair_dict = readPairs(options.pairs)

    coords_dict = readCoords(options.coords_file)

    pair_class_dict,pair_count = assessPairs(pair_dict,coords_dict,options.distance,options.error,scaffold_length_dict,\
        ref_length_dict,options.output,options.verbose)

    reportStats(pair_class_dict,pair_count,options.output)

    return 0

if __name__ == "__main__":			# this will call main and start the script
    sys.exit(main())

