#!/usr/bin/env python

#A utitlity class with a set of tools for scaling lists large numbers

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import sys
from SimpleStat import *
from ChartUtilities import ChartUtilities

#A utitlity class with a set of tools for scaling lists large numbers
class ScaleUtilities:
    """A utitlity class with a set of tools for scaling lists large numbers"""

    #uses the number of standard deviations from mean to scale a list
    #of large numbers to list of smaller numbers
    #input is a list of large numbers [x,y,z....]
    def using_deviations_from_mean(self,data):
        """uses the number of standard deviations from mean to scale a list"""
        scaled_data=[]
        stats=SimpleStat(data)
        mean = stats.get_mean()
        std = stats.get_standard_deviation()
        for num in data:
            stds_from_mean = (num-mean) / std
            scaled_data.append(stds_from_mean)
        return mean,std,scaled_data


    #takes in a list of numbers and a target range the numbers should be
    #squashed to. The returns a list of numbers that are inside the target
    #range
    #number list = list of numbers ex:  [1,33,22,22]
    #target range = a list of 2 numbers [min,max] ex: [0,100]
    #Above will squash all numbers in number list to between 0 and 100
    def scale_to_range(self,number_list,target_range):
        """squashes given list of numbers to target range and creates a new list"""
        for number in number_list:
            stats = SimpleStat(number_list)
            raw_list_max = stats.get_max()
            raw_list_min = stats.get_min()
            target_range_max=target_range[1]
            target_range_min=target_range[0]
            scaled=[]
            for raw_number in number_list:
                scaled_number = float(raw_number * ( float(target_range_max-target_range_min)  /  float(raw_list_max-raw_list_min) ))
                scaled.append(scaled_number)
        return scaled




#----------------tests------

#scale=ScaleUtilities()
#cu  = ChartUtilities()

#big_nums = [1010101,121222,121211,121223,122333]
#mean,std,scaled_nums = scale.using_deviations_from_mean(big_nums)
#print "big",big_nums
#print "scaled",scaled_nums
#print "mean",mean
#print "std",std

#method: scale_to_range()
#nums=[120,140,440,880,900,220,450,680,340]
#range=[0,100]
#scaled=scale.scale_to_range(nums,range)
#print nums
#print scaled


#all_data = [
#    [0,1,2,3,4,5,6,7,8],
#    nums,
#    [0,1,2,3,4,5,6,7,8],
#    scaled
#    ]

#cu.gen_multiple_scatters_on_same_plot(
#    all_data,
#    "index",
#    "scaled and non scaled number",
#    "Scaled Test",
#    ['Non Scaled','Scaled'],
#    'ScaleTest')




    

