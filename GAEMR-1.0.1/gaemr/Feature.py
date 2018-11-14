#
# Class for generic features (e.g. can be coverage for a window,
# gc in a window, etc.)
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

class Feature:
    """ This class represents a feature of something."""

    def __init__(self, start, end, value=None, name=None):
        self.start = int(start)
        self.end = int(end)
        self.value = value
        self.name = name

    # pretty print, can change the print
    def __str(self):
        return "%s-%s: %s (%s)" % (self.start, self.end, self.name, self.value)

    # return length of feature
    def get_length(self):
        """Return length of feature"""
        return self.end - self.start

    # start of feature
    def get_start(self):
        """Return start of feature"""
        return self.start

    # end of feature
    def get_end(self):
        """Return end of feature"""
        return self.end

    # value stored for feature
    def get_value(self):
       """Return feature's value"""
       return self.value

    # name for feature
    def get_name(self):
        """Return feature's name"""
        return self.name

    # setter functions if we want to change start and stop
    def set_start(self,new_start):
        if new_start > self.start and new_start < self.end:
            self.start = start

    def set_end(self,new_end):
        if new_end < self.end and new_end > self.start:
            self.end = new_end

    def set_name(self,new_name):
        self.name = new_name

    def set_value(self,new_value):
        self.value = new_value

    def increment_value(self,new_value):
        self.value += new_value