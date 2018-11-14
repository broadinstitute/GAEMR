#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.


import matplotlib
import matplotlib.pyplot as plt

class Chart:
    """A wrapper to reduce dependency to a certain graphing library"""
    
    #create chart
    def __init__(self,chrt,fig):
        """initialize class"""
        self.chart = chrt
        self.fig=fig
        self.curr_lib = "matplotlib"

    #save as file
    def save_as(self,name):
        """saves image to file system"""
        self.chart.savefig(name,dpi=100)
        self.chart.clf()

    #return raw plt object
    def get_raw_chart_obj(self):
        """returns chart object as is"""
        return self.chart

    #return raw figure object
    def get_raw_figure_obj(self):
        """returns fig object as is"""
        return self.fig

    #String representation of object
    def __str__(self):
        """string representation of object"""
        return "current library = "+self.curr_lib + " : object = " + str(self.chart)
        
        
    

    
