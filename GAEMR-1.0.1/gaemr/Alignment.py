#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.


import os
import subprocess
import sys
import re
from BamFile import BamFile
import PlatformConstant as pc
constant = pc.PlatformConstant()

def command(cmd):
    print cmd
    subprocess.check_call(cmd, shell=True)

class Alignment(object):
    def __init__(self,query,reference,paired=True,short=True,threads=2,output_header="aligned",output_path=".", force_unpaired=False):
         self.query = query
         self.reference = reference
         self.paired = paired
         self.short = short
         self.threads = threads
         self.output_path = output_path
     #    self.output_path = os.path.dirname(os.path.abspath(output_path))
         self.output_header = os.path.basename(output_header)
         self.force_unpaired = force_unpaired

    def is_paired(self):
         return self.paired

    def is_short(self):
         return self.short

    def is_forced_unpaired(self):
        return self.force_unpaired

    def build_SeqDict(self):
        (ref, ext) = os.path.splitext(os.path.abspath(self.reference))
        if os.path.exists(ref+".dict"):
            os.unlink(ref + ".dict")
        command("java -jar " + constant.PICARD +"/CreateSequenceDictionary.jar  R= " + self.reference +" O= " + ref +".dict TMP_DIR= " + constant.PICARD_TMP_DIR)

    def mark_duplicates(self,bam):
        output_header = bam.rstrip(".bam")
        command("java -jar " + constant.PICARD +"/MarkDuplicates.jar I= " + bam + " O= " + output_header + ".duplicates_marked.bam M= " + output_header + ".duplicates_marked.metrics  TMP_DIR= " + constant.PICARD_TMP_DIR)


    def create_bam_index(self,bam):
        command(constant.SAMTOOLS + " index " + bam)

    def merge_fastq_files(self, list):

        output = re.sub("read1","merged",list[0])
        command("cat " + list[0] + " " + list[1] + " > " + output)
        return[output]

    def convert_to_fastq(self,tmp_dir):
        if self.is_paired() or self.is_forced_unpaired():
            direction = "fr"
        else:
            direction = None

        bam = BamFile(self.query, "bam", tmp_dir + "/" + self.output_header + ".bam",direction)
        (code, fastq_list) = bam.convert_to_fastq()
        if not code:
            print "Conversion to Fastq failed"
        if self.is_forced_unpaired():
            fastq_list = self.merge_fastq_files(fastq_list)

        return fastq_list

class BwaAlignment(Alignment):
    def __init__(self,query,reference,output_header,paired=True,short=True,threads=2,forced=False):
        full_path = os.path.dirname(os.path.abspath(output_header))
        full_header = os.path.basename(output_header)
        super(BwaAlignment,self).__init__(query,reference,paired,short,threads,full_header,full_path,forced)

    def build_index(self):

        if not os.path.exists(self.reference):
            print "ERROR: Reference file: " + self.reference + " - Does not exist!!"
            sys.exit(-1)

        if not os.path.dirname(os.path.abspath(self.reference)) == os.path.abspath(self.tmp_dir):
            ref_link = os.path.abspath(self.tmp_dir)+"/"+ os.path.basename(self.reference)
            print "Reference path does not match temp_dir - link reference to tmp dir: " + os.path.dirname(os.path.abspath(self.reference)) + " -TMP: " + os.path.abspath(self.tmp_dir)
	    if os.path.exists(ref_link):
                os.unlink(ref_link)
            os.symlink(os.path.abspath(self.reference), ref_link)
            self.reference = ref_link

        if os.path.getsize(self.reference)  < 2000000000:
            algorithm = "is"
        else:
            algorithm = "bwtsw"

        command(constant.BWA + "/bwa index -a " + algorithm + " " + self.reference)

    def __run_BWA_sam(self):
        bam = self.query
        if self.is_paired():
            command(constant.BWA +"/bwa sampe -a 100000 -t " + str(self.threads) +" -f " + self.tmp_dir + "/" + self.output_header + "_bwa.sam " + self.reference + " " + self.tmp_dir + "/" + self.output_header + ".1.sai " + self.tmp_dir + "/" + self.output_header + ".2.sai " + bam + " " + bam )
        else:
            command(constant.BWA +"/bwa samse -t " + str(self.threads) + " -f " + self.tmp_dir + "/" + self.output_header + "_bwa.sam " + self.reference + " " + self.tmp_dir + "/" + self.output_header + ".0.sai " )
                
    def __run_BWA_aln(self,mate):
        # remove -q 5 from defaults. --bruce 9/14/12
        command(constant.BWA +"/bwa aln " + self.reference + " -l 32 -k 2 -t " + str(self.threads) + " -o 1 -f " + self.tmp_dir + "/" + self.output_header + "." + str(mate) + ".sai -b"+ str(mate) +" "+ self.query)

    def __add_unmapped_reads(self, build_index=False):

        if build_index:
            self.build_SeqDict()
        paired = "true"
        if not self.is_paired():
            paired = "false"
        long_read_args = ""
        if not self.is_short():
            long_read_args = " MAX_GAPS= -1 "
        command("java -jar " + constant.PICARD +"/MergeBamAlignment.jar R= " + self.reference + " UNMAPPED_BAM= " + self.query + " ALIGNED_BAM= " + self.output_path + "/" + self.output_header + "_bwa.sam OUTPUT= " + self.output_path + "/" + self.output_header + ".bam ORIENTATIONS=FR  PROGRAM_RECORD_ID=BWA PROGRAM_GROUP_VERSION=0.5.9 PROGRAM_GROUP_COMMAND_LINE=tk VALIDATION_STRINGENCY=SILENT CLIP_ADAPTERS=false ALIGNER_PROPER_PAIR_FLAGS=true PE= "+ paired + " TMP_DIR= " + self.tmp_dir + long_read_args)

        #Create BAM Index
        #self.create_bam_index(self.output_path  + "/" + self.output_header + ".bam")

    def __run_BWA_bwasw(self):
        query_fastq = self.convert_to_fastq(self.tmp_dir)
        command(constant.BWA + "/bwa bwasw -t " + str(self.threads) + " -f " + self.tmp_dir + "/" + self.output_header + "_bwa.sam " + self.reference + " " +  query_fastq[0])

        #Convert to BAM
        #aligned_sam = BamFile(self.tmp_dir + "/" + self.output_header + "_bwa.sam", "sam", self.tmp_dir + "/" + self.output_header + "_unsorted.bam")
        #aligned_sam.convert_to_bam()

        #Sort BAM
        #aligned_bam = BamFile(self.tmp_dir + "/" + self.output_header + "_unsorted.bam", "bam", self.output_path + "/" + self.output_header)
        #aligned_bam.sort_bam()

        #self.create_bam_index(self.output_path + "/" + self.output_header + ".bam")

    def clean_up_files(self):
        if self.is_short():
            if self.is_paired():
                os.unlink(self.tmp_dir + "/" + self.output_header + ".1.sai")
                os.unlink(self.tmp_dir + "/" + self.output_header + ".2.sai")
            else:
                os.unlink(self.tmp_dir + "/" + self.output_header + ".0.sai")
        else:
            if self.is_paired():
                os.unlink(self.tmp_dir + "/" + self.output_header + "_unsorted.bam")

        os.unlink(self.tmp_dir + "/" + self.output_header + "_bwa.sam")


       # os.unlink(self.tmp_dir + "/*fastq")


    def run_alignment(self, build_index= True, tmp_dir=".", mark_dups=False):
         self.tmp_dir = tmp_dir
         if build_index:
             self.build_index()
         if self.is_paired() and self.is_short():
              self.__run_BWA_aln(1)
              self.__run_BWA_aln(2)
              self.__run_BWA_sam()
              self.__add_unmapped_reads(build_index)
              final_bam = self.output_path + "/" + self.output_header + ".bam"
              if mark_dups:
                  self.mark_duplicates(final_bam)
                  final_bam = self.output_path + "/" + self.output_header + ".duplicates_marked.bam"
              self.create_bam_index(final_bam)
         elif not self.is_paired() and self.is_short():
             self.__run_BWA_aln(0)
             self.__run_BWA_sam()
             self.__add_unmapped_reads(build_index)
             final_bam = self.output_path + "/" + self.output_header + ".bam"
             if mark_dups:
                 self.mark_duplicates(final_bam)
                 final_bam = self.output_path + "/" + self.output_header + ".duplicates_marked.bam"
             self.create_bam_index(final_bam)
         elif not self.is_paired() and not self.is_short():
              self.__run_BWA_bwasw()
              self.__add_unmapped_reads(build_index)
              final_bam = self.output_path + "/" + self.output_header + ".bam"
              self.create_bam_index(final_bam)
         else:
              print "Unable to run BWA on paired: " + str(self.paired) + " and short: " + str(self.short) + " Try using BowTie"
              sys.exit(-1)



class BowTieAlignment(Alignment):
    def __init__(self,query,reference,output_header,paired=True,short=True,threads=2):
        full_path = os.path.dirname(os.path.abspath(output_header))
        full_header = os.path.basename(output_header)
        super(BowTieAlignment,self).__init__(query,reference,paired,short,threads,full_header,full_path)
        self.index_name = os.path.basename(self.reference)

    def build_index(self):
        command(constant.BOWTIE + "/bowtie2-build " + self.reference + " " + os.path.basename(self.reference))

    def build_bowtie_cmd(self, fastq_list, options=""):

        if len(fastq_list) == 1:
            read_string = "-U " + fastq_list[0]
        elif len(fastq_list) == 2:
            read_string = "-1 " + fastq_list[0] + " -2 "+ fastq_list[1]
        else:
            print "There are "+ str(len(fastq_list)) + "fastq files in the list.  There should be 1 for an unpaired BAM or 2 for a paired BAM"
            sys.exit(-1)

        command = constant.BOWTIE + "/bowtie2 -p " + str(self.threads) + " " + options + " -x " + self.index_name + " " + read_string + " -S " + self.output_path + "/" + self.output_header + "_bowtie.sam"

        return command

    def create_sorted_bam(self):
        aligned_sam = BamFile(self.tmp_dir + "/" + self.output_header + "_bowtie.sam", "sam", self.tmp_dir + "/" + self.output_header + "_unsorted.bam")
        aligned_sam.convert_to_bam()

        aligned_bam = BamFile(self.tmp_dir + "/" + self.output_header + "_unsorted.bam", "bam", self.output_path + "/" + self.output_header )
        aligned_bam.sort_bam()

        self.create_bam_index(self.output_path + "/" + self.output_header + ".bam")

    def clean_up_files(self):
        #os.unlink(self.tmp_dir + "/*fastq")
        os.unlink(self.tmp_dir + "/" + self.output_header + "_unsorted.bam")
        os.unlink(self.tmp_dir + "/" + self.output_header + "_bowtie.sam")

    def run_alignment(self,build_index=True,tmp_dir="."):
        self.tmp_dir = tmp_dir
        if build_index:
            self.build_index()
        fastq_list = self.convert_to_fastq(tmp_dir)
        if self.is_paired() and self.is_short():
            options = ""
        elif not self.is_paired() and self.is_short():
            options = ""
        elif not self.is_paired() and not self.is_short():
            options = "--local"
        elif self.is_paired() and not self.is_short():
            options = ""
        else:
            print "Unable to run BowTie on paired: " + str(self.paired) + " and short: " + str(self.short) + " Try using BWA"
            sys.exit(-1)

        command(self.build_bowtie_cmd(fastq_list, options))
        self.create_sorted_bam()
        #self.clean_up_files()




