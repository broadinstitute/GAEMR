#
# Class for manipulating sff files
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
import subprocess
from RunCommand import RunCommand
import PlatformConstant as pc

constant = pc.PlatformConstant()

class SffFile():
    """ This class represents a manipulation of SFF files. """

    def __init__(self, file, direction=None, output=None):
        self.file = file
        self.paired = 0
        if direction:
            self.paired = 1
            self.direction = direction
        if output:
            self.output = output
        else:
            self.output = self.__get_output_filename()

    def __get_output_filename(self):
        return str(self.__get_base_header()) + ".unmapped.bam"
    
    def __get_base_header(self):
        names = self.get_input_file().split('.')
        return names[0] + "." + names[1] + "." + names[2]

    def get_output_file(self):
        return self.output
    
    def convert_to_unmapped_bam(self, java_heap='32G'):
        cmd = self.__build_command(java_heap)
        if cmd:
            rc = RunCommand(cmd)
            out = rc.run_command()
        return 1

    def __build_command(self, java_heap='32G'):
        convert_command = constant.PICARD + "SffToSam.jar"
        tmp_dir = "TMP_DIR=" + constant.PICARD_TMP_DIR
        run,region = self.__get_run_region(self.get_input_file())
        return ["java","-Xmx" + java_heap, "-jar",convert_command,
                "I=",self.get_input_file(),"O=",self.get_output_file(),
                tmp_dir,"SAMPLE=",run,"RG=", self.__get_base_header()]

    def __get_run_region(self,file):
        paths = file.split("/")
        base = paths[-1].split(".")
        return paths[-1],str(1)
        
    def get_input_file(self):
        return self.file

    def is_paired(self):
        return self.paired

    def get_direction(self):
        return self.direction
