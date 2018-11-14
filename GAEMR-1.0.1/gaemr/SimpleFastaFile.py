
# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.


import os
import sys
from string import maketrans


class SimpleFastaFile:
    """ a suite of tools to manipulate simple fasta files"""

    #initialize fields
    def __init__(self,seqName):
        if not len(seqName):
            print "please enter a sequence file"
            sys.exit()
        self.sequence = dict()
        self.order = []
        self.seqFileName = seqName
        if self.exists():
            self.readSeqFile()

    #Add sequence
    def addSeq(self,name,seq):
        seq.replace("\n", "")
        if self.sequence.__contains__(name):
            self.sequence[name] = self.sequence[name] + seq.rstrip()
        else:
            self.order.append(name)
            self.sequence[name] = seq.rstrip()

    #read a sequence file
    def readSeqFile(self):
        seqFile = open(self.seqFileName,"r")
        name = ""
        for line in seqFile:
            if ">" == line[0]:
                header = line.rstrip().split(" ")
                nameField = header[0]
                name = nameField[1:]
            else:
                self.addSeq(name,line)

    #write a sequence file
    def writeSeqFile(self):
        if self.exists():
            print "File ("+self.seqFileName+") already exists.  Will not write."
            sys.exit()
        seqFile = open(self.seqFileName,"w")
        for contig in self.order:
            seqFile.write(">"+contig+"\n")
            seqFile.write(self.sequence[contig]+"\n")

    #check if file exists
    def exists(self):
        if os.path.isfile(self.seqFileName):
            return 1
        else:
            return 0

    def reverseSeq(self,name):
        seq = self.sequence[name]
        rev_seq = seq[::-1]
        trantab = maketrans("actgACTGnN","tgacTGACnN")
        rc_seq = rev_seq.translate(trantab)
        self.sequence[name] = rc_seq

    #check if unique
    def isUnique(self,query):
        for name,seq in self.sequence.iteritems():
            if query in seq:
                return 0
        return 1

    #extract right sequence tag
    def extractSeqFromRightExtreme(self,seq,length):
        if len(seq)>length:
            start = len(seq)-length
        else:
            start = 0
        return seq[start:]

    #extract left sequence tag
    def extractSeqFromLeftExtreme(self,seq,length):
        if length < len(seq):
            end = length
        else:
            return seq[0:]
        return seq[:end]

    #extract sequence tag
    def extractSequenceTags(self):
        seqTags =[]
        for fasta in self.sequence.values():
            seqTags.append(self.extractSeqFromLeftExtreme(fasta, 20))
            seqTags.append(self.extractSeqFromRightExtreme(fasta, 20))
        return seqTags

    #extract kmers
    def extractKmers(self,kMerSize=40):
        return extractKmersToFile(None,kMerSize,False)

    #extract kmers
    def extractKmersToFile(self,out_file_name=None,kMerSize=40,write=True):
        kMer_list = []
        kMers = sys.stdout
        if out_file_name:
            kMers = open(out_file_name,"w")
        for contig_name,fasta in self.sequence.iteritems():
            for i in range(0,len(fasta) - kMerSize + 1):
                if write:
                    header = (">"+contig_name.replace("\n","") + "_" + str(i) + "\n").replace("|","_").replace(" ","__")
                    kMers.write(header)
                    kMers.write(fasta[i:i+kMerSize].upper() + "\n")
                else:
                    kMer_list.append(fasta[i:i+kMerSize].upper())
        kMers.close()
        if write:
            return out_file_name
        return kMer_list


    #extract kmers to dictionary, key is header, value is fasta
    def get_kmers_structure(self,kmer_size=40,write=True):
        kmer_struct =dict() 
        for contig_name,fasta in self.sequence.iteritems():
            for i in range(0,len(fasta) - kmer_size + 1):
                if write:
                    header = (">"+contig_name.replace("\n","") + "_" + str(i) + "\n").replace("|","_").replace(" ","__")
                    kmer_struct[header] =fasta[i:i+kmer_size].upper()
        return kmer_struct




    #public method gets access to sequence
    def get_sequence(self):
        if not len(self.sequence):
            raise "ERROR: sequence not populated"
        return self.sequence

