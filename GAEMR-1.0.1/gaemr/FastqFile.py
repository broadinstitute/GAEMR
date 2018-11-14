#
# Class for manipulating fastq files
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
from RunCommand import RunCommand
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import PlatformConstant as pc

constant = pc.PlatformConstant()

class FastqFile():
    """ This class represents a manipulation of Fastq files. """

    def __init__(self, files, sample=None, qual_format=None, direction=None, output=None):
        self.files = []
        self.sample = 'sample'
        if sample:
            self.sample = sample

        if qual_format:
            self.qual_format = qual_format

        for i in files:
            self.files.append(i)
        if qual_format:
            self.qual_format = qual_format
        else:
            self.qual_format = "Standard"

        if output:
            self.output = output
        else:
            self.output = self.__get_output_filename(files[0])
        if direction:
            self.direction = direction.upper()
            self.__check_direction()
        
    def __get_output_filename(self,name):
        return str(self.__get_base_header(name)) + ".unmapped.bam"
    
    def __get_base_header(self,name):
        names = name.split('.')
        return names[0]

    def convert_to_unmapped_bam(self, java_heap='32G'):
        cmd = self.__build_command(java_heap)
        if cmd:
            rc = RunCommand(cmd)
            out = rc.run_command()
        return 1
    
    def __build_command(self,java_heap):
        convert_command = constant.PICARD + "FastqToSam.jar"
        tmp_dir = "TMP_DIR=" + constant.PICARD_TMP_DIR
        sample_name = "SM=" + self.get_sample_name()
        qual_arg = "V=" + self.get_qual_format()

        num_files = self.__get_num_files()
        if num_files == 1:
            return ["java","-Xmx" + java_heap, "-jar",convert_command, "F1=",self.__get_single_file(0), "O=", self.get_output_file(),
                    sample_name,qual_arg,tmp_dir]
        if num_files == 2:
            return ["java","-Xmx" + java_heap, "-jar",convert_command, "F1=",self.__get_single_file(0), "F2=",self.__get_single_file(1),"O=",
                    self.get_output_file(),sample_name,qual_arg,tmp_dir]
        return 0

    def __get_file_inputs(self,index,file):
        return "F" + str(index + 1) + "=" + file
            
    def __get_single_file(self,index):
        if index < self.__get_num_files():
            return self.files[index]
        return None
    
    def get_files(self):
        return self.files

    def __check_direction(self):
        if self.__get_num_files() > 1:
            direction = self.get_pairing()
            if direction:
                direction = direction.upper()

                # reverse both
                if direction == 'RF':
                    for i in xrange(int(self.__get_num_files())):
                        file = self.__get_single_file(i)
                        out = self.__get_rc_fastq_name(file)
                        self.reverse_file(file,out)
                        self.__set_file(i,out)

                # reverse 2nd
                if direction == 'FF':
                    infile = self.__get_single_file(1)
                    out = self.__get_rc_fastq_name(infile)
                    self.reverse_file(infile,out)
                    self.__set_file(1,out)
                    
                # reverse 1st
                if direction == 'RR':
                    infile = self.__get_single_file(0)
                    out = self.__get_rc_fastq_name(infile)
                    self.reverse_file(infile,out)
                    self.__set_file(0,out)
                    
    def __get_rc_fastq_name(self,file):
        return re.sub("fastq","rc.fastq",file)
    
    def reverse_file(self,infile,outfile):
        input = open(infile, 'r')
        output = open(outfile, 'w')
        for record in SeqIO.parse(input,'fastq'):
            s = record.seq
            q = record.letter_annotations['phred_quality']
            del record.letter_annotations['phred_quality']
            record.seq = s.reverse_complement()
            q.reverse()
            record.letter_annotations['phred_quality'] = q
            SeqIO.write(record,output,'fastq')

    def get_output_file(self):
        return self.output

    def get_pairing(self):
        return self.direction

    def __set_file(self,index,name):
        self.files[index] = name
        
    def __get_num_files(self):
        return len(self.files)

    def get_qual_format(self):
        return self.qual_format
    
    def get_sample_name(self):
        return self.sample
