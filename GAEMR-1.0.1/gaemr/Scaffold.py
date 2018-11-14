#
# Class for holding scaffold information from a SeqIO object
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

from Bio import SeqIO
from Feature import Feature 

class Scaffold:
    """ This class represents a scaffold and its constituent parts."""

    def __init__(self, scaffold, minGapSize, minConSize, contig_ids=None):
        self.scaffold = scaffold
        self.minGapSize = minGapSize
        self.minConSize = minConSize
        self.contigs = []
        self.seq_start = 0
        self.seq_end = len(scaffold)
        self.__scaffold_contigs(contig_ids)

    # private function that gets the contigs based on Ns and gap sizes
    def __scaffold_contigs(self, contig_ids=None):
        seq = str(self.get_seq()).upper()
        s_id = self.get_name()
        slen = len(seq)

        i = c_start = 0

        contig_count = 0
        value = 1
        id = None
        last_contig = 0

        while True:
            i = seq.find('N', i)
            if i < 0: break
            # count consecutive Ns
            n_start = i
            while i < slen and seq[i] == 'N':
                i += 1
            # this many Ns in a row constitute a contig break (gap)
            n_len = i - n_start
            if n_len >= self.minGapSize:
                c_len = n_start - c_start
                if c_len >= self.minConSize:
                    id = s_id + "_c" + str(contig_count + 1)
                    if contig_ids:
                        id = contig_ids[contig_count]
                    self.contigs.append(Feature(c_start,n_start,value,id))
                    contig_count += 1
                    last_contig = n_start
                elif contig_count == 0:
                    self.seq_start = i
                c_start = i
                #contig_count += 1

        if last_contig < slen:
            if slen - c_start > self.minConSize:
                id = s_id + "_c" + str(contig_count + 1)
                if contig_ids:
                    id = contig_ids[contig_count]
                self.contigs.append(Feature(c_start,slen,value,id))
            else:
                self.seq_end = last_contig

        self.get_contig_lengths_list()
        assert self.get_contig_length() + self.get_gap_length() == self.get_length()

    # public funtion to get the ungapped length of scaffold 
    def get_contig_length(self):
        """Returns the total contig length for the scaffold."""
        sum = 0
        for f in self.contigs:
            sum += f.get_length()
        return sum

    # public function to get contig sizes
    def get_contig_lengths_list(self):
        """Returns the contig sizes in a list."""
        contig_sizes = []
        for f in self.contigs:
            contig_sizes.append(f.get_length())
        return contig_sizes

    # public function to get the non-contig sequence
    def get_gaps_lengths_list(self):
        """Returns a list of non-contig (gaps) features."""
        end = -1
        gaps = []
        for f in self.contigs:
            if end == -1:
                start = f.get_start()
                end = f.get_end()
                #if start:
                #    gaps.append(f.get_start())
            else:
                gaps.append(f.get_start() - end)
                end = f.get_end()

        return gaps

    # public function to get total gap length
    def get_gap_length(self):
        """Returns sum of gap lengths."""
        sum = 0
        gaps = self.get_gaps_lengths_list()
        for f in gaps:
            sum += f
        return sum

    # public function to get the total GC of the super
    def get_GC(self):
        """Returns total gc in contigs in scaffolds"""
        gc = total = 0
        for f in self.contigs:
            seq = self.scaffold.seq[f.get_start():f.get_end()].upper()
            gc += seq.count('G') + seq.count('C')
            total += f.get_length()
        return gc,total

    # public function to get length of scaffold
    def get_length(self):
        """Returns scaffold's length."""
        return len(self.scaffold.seq[self.seq_start:self.seq_end])

    # public function to get name of scaffold
    def get_name(self):
        """Returns scaffold's id"""
        return self.scaffold.id

    # public function to set name of scaffold
    def set_name(self, name):
        """Sets scaffold's id"""
        self.scaffold.id = name

    # public function to get sequence of scaffold
    def get_seq(self):
        """Returns scaffold's sequence."""
        return self.scaffold.seq[self.seq_start:self.seq_end]

    def trim_seq(self,trim_start,trim_stop):
    #Calculate portion of scaffold to keep and then trim
        start = 0
        end = 0
        if not trim_start:
            start = trim_stop + 1
            end = self.get_length()
        if trim_stop == self.get_length():
            start = 0
            end = trim_start - 1

        trimmed_seq = self.get_seq()[start:end]
        self.scaffold.seq = trimmed_seq

    def get_contig_sequences(self):
        contigs = {}
        contigs_list = []
        for f in self.contigs:
            contigs_list.append(f.get_name())
            contigs[f.get_name()] = self.scaffold.seq[f.get_start():f.get_end()].upper()
        return contigs,contigs_list

    def get_contig_ids(self):
        contig_ids = []
        for f in self.contigs:
            contig_ids.append(f.get_name())
        return contig_ids