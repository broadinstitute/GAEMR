#
# Class for manipulating fasta files
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
import os
from RunCommand import RunCommand
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from FastqFile import FastqFile
import PlatformConstant as pc

constant = pc.PlatformConstant()

class FastaFile():
    """ This class represents a manipulation of Fasta files. """

    def __init__(self, files, sample, default_qual, quals=None, direction=None, output=None):
        self.files = []
        self.quals = []
        self.fastqs = []
        self.sample = sample
        self.qual_format = "Standard"
        self.default_qual = default_qual
        
        for i in files:
            self.files.append(i)
        if quals:
            for i in quals:
                self.quals.append(i)

        self.direction = 'fr'
        if direction:
            self.direction = direction
            
        if output:
            self.output = output
        else:
            self.output = self.__get_output_filename(files[0])

        self.fastqs = self.__convert_to_fastq(self.files,self.quals)
        self.fq_fileobject = FastqFile(self.fastqs,self.sample,self.qual_format,self.direction,self.output)
        
    def __get_output_filename(self,name):
        return str(self.__get_base_header(name)) + ".unmapped.bam"
    
    def __get_base_header(self,name):
        names = name.split('.')
        return names[0]

    def __get_fastq_name(self,file):
        return re.sub("fasta","fastq",file)

    def __convert_to_fastq(self,files,quals):
        fastqs = []

        for i in range(int(self.__get_num_files())):
            fastq_file = self.__get_fastq_name(files[i])
            out_handle = None
            if quals:
                rec_iter = PairedFastaQualIterator(open(self.__get_single_file(i,files), "rU"),
                                                   open(self.__get_single_file(i,quals), "rU"))
                out_handle = open(fastq_file, "w")
                SeqIO.write(rec_iter, out_handle, "fastq")

            else:
                f_handle = open(self.__get_single_file(i,files),"rU")
                out_handle = open(fastq_file, "w")
                for record in SeqIO.parse(f_handle,"fasta"):
                    qual_string = self.__get_default_qual() * len(record.seq)
                    out_handle.write('@%s\n%s\n+%s\n%s\n' % (record.id, record.seq,"",qual_string))
                    
            fastqs.append(fastq_file)
            out_handle.close()
        return fastqs

    def __get_default_qual(self):
        if int(self.default_qual) > constant.QUAL_MAX:
            return chr(constant.QUAL_MAX)
        return chr(int(self.default_qual) + 33)
 
    def convert_to_unmapped_bam(self, java_heap='32G'):
        return_value = self.fq_fileobject.convert_to_unmapped_bam(java_heap)
        self.__remove_tmp_fastq_files()
        return return_value

    def __get_fastqs(self):
        return self.fastqs

    def __remove_tmp_fastq_files(self):
        for i in self.__get_fastqs():
            os.remove(i)

    def __get_single_file(self,index,files):
        if index < self.__get_num_files:
            return self.files[index]
        return None

    def __set_file(self,index,name):
        self.files[index] = name
        
    def __get_num_files(self):
        return len(self.files)

    def get_files(self):
        return self.files
