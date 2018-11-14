#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

from numpy import array

#a simple stat class around numpy
class SimpleStat():
    """A wrapper class to simplify some basic statistical analysis"""

    #init class: convert input list into numpy array
    def __init__(self,my_list):
        """start class: covert input list of numbers into numpy array"""
        self.data_set = array((my_list))

    #calculate mean
    def get_mean(self):
        """get mean"""
        mean=0.0
        try:
            mean = self.data_set.mean()
        except:
            raise "cannot calculate mean from data"
        return mean

    #standard deviation
    def get_standard_deviation(self):
        """get standard deviation"""
        std=0.0
        try:
            std = self.data_set.std()
        except:
            raise "cannot calculate standard deviation from data"
        return std

    #get biggest number
    def get_max(self):
        """find the largest number in list"""
        return self.data_set.max()

    #get smallest number
    def get_min(self):
        """find the smallest number in list"""
        return self.data_set.min()
