#
# Class to handle running commands and checking errors
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

class RunCommand:
    """ This class represents running a command and checking errors. """

    def __init__(self,command):
        self.command = command

    # public function to run the command and check for errors
    def run_command(self,check_output=1,check_errors=1):
        """Runs a command and gets errors if any."""
        command = self.command
        p = subprocess.Popen(command, stdout=subprocess.PIPE)
        if check_output:
            out, err = p.communicate()
            if p.returncode and check_errors:
                raise subprocess.CalledProcessError(p.returncode, command)
            return out
        else:
            return p

    # private function to make the command print out pretty
    def __to_string(self):
        return ''.join([str(i) + " " for i in self.command])

    # public function to get the command to run
    def get_command(self):
        """Returns command string"""
        return self.__to_string()
