#!/usr/bin/env python

import sys
import re
import os
from gaemr.RunCommand import RunCommand
from optparse import OptionParser
from gaemr.SimpleFastaFile import SimpleFastaFile
from gaemr.CoordsFile import CoordsFile
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options] nucmer.delta")

parser.add_option('--output', '-o', action="store", default="tiling_ono", type='string', dest="output",
                  help='Output prefix.')
parser.add_option('--nounaligned', '-n', action="store_true", default=False, dest="nounaligned",
                  help='Do not print unaligned contigs at the end of the file.')
parser.add_option('--circular', '-c', action="store_true", default=False, dest="circular",
                  help='Assume reference sequences are circular, and cut at reference start when spanning.')
parser.add_option('--rename', '-r', action="store_true", default=False, dest="rename",
                  help='Rename ono\'d contigs and scaffolds (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply a delta file from a NUCMER alignment.")
    sys.exit(-1)

def __build_nucmer_command(ref,seq,arg_list):
    return [constant.NUCMER] + arg_list + [ref] + [seq]

def __build_showtiling_command(delta,arg_list):
    return [constant.SHOWTILING] + arg_list + [delta]

def main():

    delta = args[0]
    delta_file = open(delta,'r')
    reference,query = delta_file.readline().rstrip().split(" ")
    query_fasta = SimpleFastaFile(query)
    query_seqs = query_fasta.get_sequence()

    if os.path.isfile(options.output+".interim.fasta"):
        remove_command = "rm",options.output+".interim.fasta"
        rc = RunCommand(remove_command)
        print rc.get_command()
        rc.run_command()

    ono_fasta = SimpleFastaFile(options.output+".interim.fasta")

    print "Ordering and orienting using", delta
    print "Reference", reference
    print "Query", query

    arg_list = ['-R', '-v 5', '-V 0', '-u', options.output+'.unplaced']
    if options.circular:
        arg_list += ["-c"]
    rc = RunCommand(__build_showtiling_command(delta,arg_list))
    print "Executing:",rc.get_command()
    print ""
    
    left_cut = 0
    
    for line in rc.run_command().splitlines():
    
        wraps = False
        if (line.startswith(">")):

            #IF THERE WAS A WRAP ON THE PREVIOUS REFERENCE SEQUENCE
            if (left_cut):
                print left_text
                ono_fasta.addSeq(left_query+"_left",left_side)
                ono_fasta.reverseSeq(left_query+"_left")
                query_seqs[left_query] = 0
                
            ref,bases,null = line.rstrip().split(" ")
            print "Ordering to", ref
            left_cut = 0

        else:
            
            rstart,rend,qstart,qend,cov,id,ori,query = line.rstrip().split("\t")
            #CHECK TO SEE IF THIS IS AN ALIGNMENT WHICH WRAPS AROUND THE END OF A REFERENCE SEQUENCE    
            if (int(rstart) < 0 and options.circular):
                #print "WRAP", line.rstrip().split("\t")
                wraps = True
                
                arg_list = ['-R', '-a', '-v 5', '-g -1', '-V 0', '-u', options.output+'.unplaced','-c']
                rc = RunCommand(__build_showtiling_command(delta,arg_list))
                #print "Executing:",rc.get_command()
                
                for line in (rc.run_command().splitlines()):
                    align = line.rstrip().split("\t")
                    #print align
                    #print align[0], align[12], query
                    if (int(align[0]) == 1 and align[12] == query):
                        #print "CLEAN WRAP FOUND!", align
                        #PRINT RIGHT SIDE OF QUERY AND STORE LEFT CUT SITE
                        if (ori == "+"):
                            print "\tStoring left side of the query to print at the end of alignments to this reference sequence"
                            left_cut = int(align[2])-1
                            left_side = query_seqs[query][:left_cut]
                            left_query = query
                            left_text = "\tWriting left side of " + query + " (1 to " + str(left_cut) + ") aligned " + cov + " at " + id + "% identity in the forward orientation."

                            right_cut = int(align[2])-1
                            print "\tWriting right side of " + query + " (" + str(right_cut+1) + " to end) aligned " + cov + " at " + id + "% identity starting at 1 in the forward orientation."
                            right_side = query_seqs[query][right_cut:]
                            ono_fasta.addSeq(query+"_right",right_side)
                        else:
                            print "\tStoring left side of the query to print at the end of alignments to this reference sequence"
                            left_cut = int(align[2])
                            left_side = query_seqs[query][left_cut:]
                            left_query = query
                            left_text = "\tWriting left side of " + query + "(" + str(left_cut+1) +" to end) aligned " + cov + " at " + id + "% identity in the reverse orientation."
                    
                            right_cut = int(align[2])
                            print "\tWriting right side of " + query + " (1 to " + str(right_cut) + ") aligned " + cov, "at", id + "% identity starting at 1 in the reverse orientation."
                            right_side = query_seqs[query][:right_cut]
                            ono_fasta.addSeq(query+"_right",right_side)
                            ono_fasta.reverseSeq(query+"_right")
                            
                if not left_cut:
                    print "\tIt appears this query sequence overlaps, but no clean cut site was found.  Printing as normal."
                    wraps = False

            #IF THERE IS NO WRAP
            if (not wraps):
                if ori == "+":
                    print "\tWriting", query, "aligned", cov, "at", id + "% identity starting at", rstart, "in the forward orientation."
                    ono_fasta.addSeq(query,query_seqs[query])
                else:
                    print "\tWriting", query, "aligned", cov, "at", id + "% identity starting at", rstart, "in the reverse orientation."
                    ono_fasta.addSeq(query,query_seqs[query])
                    ono_fasta.reverseSeq(query)    
                query_seqs[query] = 0
                
    #IF THE LAST REFERENCE HAD A WRAP, PRINT LEFT SIDE
    if (left_cut):
        print left_text
        ono_fasta.addSeq(left_query+"_left",left_side)
        if "reverse" in left_text:
            ono_fasta.reverseSeq(left_query+"_left")
        query_seqs[left_query] = 0

    #PRINT ANY UNALIGNED SEQUENCES
    for contig in query_seqs:
        if query_seqs[contig]:
            print "Writing unaligned sequence:", contig
            ono_fasta.addSeq(contig,query_seqs[contig])

    ono_fasta.writeSeqFile()
    make_assembly_command = ["make_standard_assembly_files.py","-S",options.output+".interim.fasta","-o",options.output+".ono"]
    if options.rename:
        make_assembly_command += ['-r']
    rc = RunCommand(make_assembly_command)
    print "Executing",rc.get_command()
    rc.run_command()

if __name__ == "__main__":
    sys.exit(main())
