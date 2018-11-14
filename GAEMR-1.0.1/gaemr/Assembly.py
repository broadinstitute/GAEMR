#
# Class for manipulating assembly fasta data
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
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Scaffold import Scaffold
from AgpFile import AgpFile
import PlatformConstant as pc
constant = pc.PlatformConstant()

class Assembly(object):
    """ This class represents an assembly and its constituent parts."""

    def __init__(self, fasta=None, name=None, assembler=None, minGapSize=None, minConSize=None, minScaffSize=None, agp=None):
        self.fasta = fasta

        self.minGapSize = constant.MIN_GAP_SIZE
        self.minConSize = constant.MIN_CONTIG_SIZE
        self.minScaffSize = constant.MIN_SCAFFOLD_SIZE

        if name:
            self.name = name
        if assembler:
            self.assembler=assembler
        if minGapSize:
            self.minGapSize = minGapSize
        if minConSize:
            self.minConSize = minConSize
        if minScaffSize:
            self.minScaffSize = minScaffSize
        self.scaffolds = []

        # can set up assembly with just fasta or agp now
        if fasta:
            if agp:
                f_index = SeqIO.index(fasta,"fasta")
                self.agp = AgpFile(agp)
                self.__read_agp(f_index)
            else:
                self.__read_fasta(fasta)

            self.__clean_scaffolds()

    def __clean_scaffolds(self):
        scaffolds = []
        for s in self.scaffolds:
            if s.get_length() > self.minScaffSize:
                scaffolds.append(s)
        self.scaffolds = scaffolds

    # private function to read fasta file
    def __read_fasta(self,fasta):
        """Makes scaffolds from fasta file"""
        for scaffold in SeqIO.parse(open(fasta, 'r'), "fasta"):
            slen = len(scaffold)
            if slen < self.minScaffSize:
                continue
            self.scaffolds.append(Scaffold(scaffold,self.minGapSize,self.minConSize))

    # private function to make scaffold from agp
    def __read_agp(self,f_index):
        """Gets agp information to make a scaffold"""
        agp = self.__get_agp_object()
        scaffolds = agp.get_agp_scaffolds()
        for i in scaffolds:
            f_dict = agp.get_agp_file_record(i)
            self.__build_scaffold(agp,i,f_dict,f_index)

    # private function to build the scaffold from agp
    def __build_scaffold(self,agp,scaffold_id,f_dict,f_index):
        """Makes the scaffold records based on agp input"""
        scaffold_seq = ""
        contig_ids = []
        for i in f_dict:
            if agp.is_gap(scaffold_id,i):
                length = agp.get_feature_length(scaffold_id,i)
                gap = 'N' * int(length)
                scaffold_seq += gap
            else:
                contig_id = agp.get_contig_id(scaffold_id,i)
                sequence = f_index[contig_id].seq
                if agp.is_contig_reverse_complemented(contig_id):
                    tmp = f_index[contig_id].reverse_complement()
                    sequence = tmp.seq
                scaffold_seq += sequence
                contig_ids.append(contig_id)
        scaff_record=SeqRecord(scaffold_seq,scaffold_id)
        if len(scaffold_seq) > self.minScaffSize:
            self.scaffolds.append(Scaffold(scaff_record,self.minGapSize,self.minConSize,contig_ids))

    # private function to get the agp object
    def __get_agp_object(self):
        """Returns agp object"""
        return self.agp

    # public function to set the assembler post-creation
    def set_assembler(self, assembler):
        """Sets the assembler as some value."""
        self.assembler=assembler

    # public function to set the name of the assembly post-creation
    def set_name(self,name):
        """Sets the assembly name."""
        self.name = name

    # public function to set the minimum gap size post-creation
    def set_min_gap(self,size):
        """Sets the minimum gap size"""
        self.minGapSize = size

    # public function to set the minimum contig size post-creation
    def set_min_con_size(self,size):
        """Sets the minimum contig size"""
        self.minConSize = size

    # public function to set minimum scaffold size post-creation
    def set_min_scaff_size(self,size):
        """Sets the minimum scaffold size."""
        self.minScaffSize = size

    # public function to get assembler name
    def get_assembler(self):
        """Returns assembler name"""
        return self.assembler

    # public function to get assembly name
    def get_assembly_name(self):
        """Returns assembly name."""
        return self.name

    # public function to get the total length
    def get_gapped_length(self):
        """Returns total length of assembly in scaffolds."""
        sum = 0
        for s in self.scaffolds:
            sum += s.get_length()
        return sum

    # public function to get the length of contigs
    def get_ungapped_length(self):
        """Returns total length of contigs in assembly."""
        sum = 0
        for s in self.scaffolds:
            sum += s.get_contig_length()
        return sum

    # public function to get the total gc
    def get_assembly_GC(self):
        """Returns total gc in contigs in assembly."""
        gc = total = 0
        for s in self.scaffolds:
            (g,t) = s.get_GC()
            gc += g
            total += t
        return round((float(gc)/float(total))*100, 2)

    # private function to sort scaffolds by size
    def __get_rsorted_scaffs(self):
        """Returns reversed list of scaffold sizes"""
        sort_scaff = []
        for s in self.scaffolds:
            sort_scaff.append(s.get_length())
        sort_scaff.sort(reverse=True)
        return sort_scaff

    # private function to sort contigs by size
    def __get_rsorted_contigs(self):
        """Returns reversed list of contig sizes"""
        sort_contigs = []
        for s in self.scaffolds:
            sort_contigs += s.get_contig_lengths_list()
        sort_contigs.sort(reverse=True)
        return sort_contigs

    # private function to sort gaps by size
    def __get_rsorted_gaps(self):
        """Returns reversed list of gap sizes"""
        sort_gaps = []
        for s in self.scaffolds:
            sort_gaps += s.get_gaps_lengths_list()
        sort_gaps.sort(reverse=True)
        return sort_gaps

    # public function to get max scaffold in assembly
    def get_max_scaffold(self):
        """Return largest scaffold size"""
        return self.__get_rsorted_scaffs()[0]

    # public function to get the max contig in assembly.
    def get_max_contig(self):
        """Return largest contig size."""
        return self.__get_rsorted_contigs()[0]

    # public function to get cumulative scaffold sizes
    def get_cumulative_scaffolds(self):
        """Returns list of cumulative sizes in scaffolds."""
        return self.__cumulative_data(self.__get_rsorted_scaffs())

    # public function to get cumulative contig sizes in assembly
    def get_cumulative_contigs(self):
        """Returns list of cumulative sizes in contigs."""
        return self.__cumulative_data(self.__get_rsorted_contigs())

    # private function to calculate cumulative data
    def __cumulative_data(self,toTotal):
        """Returns cumulative value list"""
        total = 0
        cumulative = []
        for v in toTotal:
            total += v
            cumulative.append(int(total))
        return cumulative

    # private function to calculate the N* of assembly
    def __nfrac(self, lengths, fraction):
        """Returns N-something of the assembly."""
        total = sum(lengths)
        partial = count = 0
        target = total * fraction
        for length in lengths:
            partial += length
            count += 1
            if partial >= target:
                return length, count
        return 0, count

    # public function to get the N* of assembly
    def nXX(self, things, XX):
        """Returns N-something of the contigs or scaffolds."""
        length, count = self.__nfrac(things, XX / 100.0)
        return length

    # public function to get the N50 of something
    def n50(self, things):
        """Returns N50 of contigs or scaffolds"""
        return self.nXX(things, 50.0)

    # public function to get the N90 of something
    def n90(self, things):
        """Returns N90 of contigs or scaffolds."""
        return self.nXX(things, 90.0)

    # public function to get N* for scaffolds
    def scaffoldN90(self):
        return self.n90(self.__get_rsorted_scaffs())

    def scaffoldN50(self):
        return self.n50(self.__get_rsorted_scaffs())

    # public functions to get N* for contigs
    def contigN90(self):
        return self.n90(self.__get_rsorted_contigs())

    def contigN50(self):
        return self.n50(self.__get_rsorted_contigs())

    # private function to calculate mean on list
    def __list_mean(self, list):
        """Returns mean value from list"""
        total = 0
        for v in list:
            total += v
        return int(round(float(total)/float(len(list)),0))

    # public functions to get means for scaffolds and contigs
    def scaffold_mean(self):
        return self.__list_mean(self.__get_rsorted_scaffs())

    def contig_mean(self):
        return self.__list_mean(self.__get_rsorted_contigs())

    # public function to get number of scaffolds
    def num_scaffolds(self):
        """Returns number of scaffolds"""
        return len(self.get_scaffolds())

    # public function to get number of contigs
    def num_contigs(self):
        """Returns number of contigs"""
        contigs = 0
        for s in self.scaffolds:
            #print "SCAFFOLD: ",
            #print s.get_name(),
            #print s.get_length()
            #print s.get_contig_lengths_list()
            contigs += len(s.get_contig_lengths_list())
        return contigs

    # private function to see if we have gaps
    def __check_num_gaps(self):
        """Return whether or not we have gaps."""
        return self.num_gaps() > 0

    # public function to get gap number
    def num_gaps(self):
        """Returns number of gaps"""
        return self.num_contigs() - self.num_scaffolds()

    # public function to get maximum gap size
    def gap_max(self):
        """Returns maximum gap size."""
        if self.__check_num_gaps():
            return self.__get_rsorted_gaps()[0]
        return 0

    # public function to get mean gap size
    def gap_mean(self):
        """Returns gap mean size."""
        if self.__check_num_gaps():
            return self.__list_mean(self.__get_rsorted_gaps())
        return 0

    # public function to get gap N50
    def gapN50(self):
        """Return gap N50"""
        if self.__check_num_gaps():
            return self.n50(self.__get_rsorted_gaps())
        return 0

    # public function to get gap length
    def gap_length(self):
        """Return total gap length."""
        if self.__check_num_gaps():
            return self.get_gapped_length() - self.get_ungapped_length()
        return 0

    # public function to add a scaffold anew to assembly
    def add_scaffold(self,scaffold):
        """Adds scaffold to an assembly"""
        if len(scaffold.seq) > self.minScaffSize:
            self.scaffolds.append(Scaffold(scaffold,self.minGapSize,self.minConSize))
            return 1
        return 0

    # public function to delete a scaffold from assembly
    def delete_scaffold(self,scaffold_id):
        """Deletes a scaffold from assembly"""
        for s in self.scaffolds:
            name = s.get_name()
            if re.match(scaffold_id,name):
                self.scaffolds.remove(s)
                return 1
        return 0

    # public function to see if scaffold exists
    def contains_scaffold(self,scaffold_id):
        """Returns whether or not a scaffold is in assembly"""
        for s in self.scaffolds:
            name = s.get_name()
            if re.match(scaffold_id,name):
                return 1
        return 0

    # public function to break a scaffold at a position
    def break_scaffold_at_position(self,scaffold_id,position,in_place=False):
        """Breaks scaffold and puts in back in place or in a new scaffold"""
        for s in self.scaffolds:
            name = s.get_name()
            s_seq = s.get_seq()
            s_len = s.get_length()
            if re.match(scaffold_id,name):
                first_seq = SeqRecord(s_seq[:position],id=name + "_1",
                                      description="Derived from " + name + " position 1 to " + str(position))
                second_seq = SeqRecord(s_seq[position:],id=name + "_2",
                                       description="Derived from " + name +
                                       " position " +  str(position + 1) + " to " + str(s_len))
                self.delete_scaffold(scaffold_id)
                if in_place:
                    self.add_scaffold(first_seq + second_seq)
                    return 1
                self.add_scaffold(first_seq)
                self.add_scaffold(second_seq)
                return 1
        return 0

    # public function to get the scaffolds in assembly
    def get_scaffolds(self):
        return self.scaffolds

    # private function to print out fasta
    def __to_fasta(self, seq_dict, seq_list, outfile, rename, num, type="contigs"):
        for i in seq_list:
        #for i in sorted(seq_dict.iterkeys()):
            if rename:
                if type == 'contigs':
                    outfile.write(">contig%06d\n" % (num))
                else:
                    outfile.write(">scaffold%05d\n" % (num))
                num += 1
            else:
                outfile.write(">" + str(i) + "\n")
            for x in xrange(0,len(seq_dict[i]),60):
                outfile.write(str(seq_dict[i][x:x+60]) + "\n")
        return num

    # public function to print contigs from assembly
    def print_assembly_contigs(self,output=None,rename=False):
        outfile = sys.stdout
        if output:
            try:
                outfile = open(output,'w')
            except IOError as (errno,strerror):
                print "I/O Error({0}): {1}".format(errno,strerror)
                return -1
        nc = 1
        for s in self.get_scaffolds():
            c_dict,c_list = s.get_contig_sequences()
            nc = self.__to_fasta(c_dict,c_list,outfile,rename,nc,type = "contigs")
        outfile.close()

    # private helper function to make an agp from assembly - in progress
    def __to_agp(self,want_version_2=True,rename=True):
        a = AgpFile()
        a.make_agp_from_assembly(self.get_scaffolds(),want_version_2,rename)
        return a

    def print_assembly_agp(self,output=None,want_version_2=True,rename=True):
        a = self.__to_agp(want_version_2,rename)
        a.print_agp(output)

    # private helper function to get the scaffolds into the right format
    # for printing out fasta output
    def __to_dict(self,scaffolds):
        scaffold_dict = {}
        scaffold_list = []
        for i in scaffolds:
            scaffold_list.append(i.get_name())
            scaffold_dict[i.get_name()] = i.get_seq()
        return scaffold_dict,scaffold_list

    # public function to print out assembly
    def print_assembly(self,file=None,rename=False):
        """Prints out our assembly."""
        ns = 1
        s_dict,s_list = self.__to_dict(self.get_scaffolds())

        if s_dict and s_list:
            output = sys.stdout
            if file:
                output = file
            try:
                assembly = open(output,'w')
            except IOError as (errno,strerror):
                print "I/O Error({0}): {1}".format(errno,strerror)
                return -1

            ns = self.__to_fasta(s_dict,s_list,assembly,rename,ns,type = "scaffolds")
            assembly.close()
            return 0
        return 1

    # public function to print out assembly
    def write_assembly_to_file(self,file):
        """Write out our assembly to file."""
        self.print_assembly(file)
        return 0

    # public function to return scaffold by id
    def get_scaffold_by_name(self,scaffold_id):
        """Returns scaffold object from  assembly"""
        for s in self.scaffolds:
            name = s.get_name()
            if re.match(scaffold_id,name):
                return s
        return 0
