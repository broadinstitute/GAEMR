#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from numpy.random import randn
from matplotlib import cm

from gaemr.Chart import Chart
from gaemr.SimpleStat import SimpleStat
from mpl_toolkits.mplot3d import Axes3D
from gaemr.RadarChart import *
import numpy as np
from matplotlib import font_manager as fm
import matplotlib.image as image
from gaemr.PlatformConstant import *
from matplotlib.ticker import AutoMinorLocator

from matplotlib.font_manager import FontProperties

import matplotlib.gridspec as gridspec


class ChartUtilities:
    """ A tool kit for generating different types of charts """

    #generates a bar plot, data=list of numbers. Please specify ouput
    #file name without any extensions.
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_bar_chart(self,data,x_name,y_name,main_name,output_file_name,ret_obj=False):
        """generate a simple bar plot of a list of numbers"""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_title(main_name)
        ax.grid(True)
        out_name = output_file_name + ".png"
        plt.plot(data)
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:           
            plt.savefig(out_name,dpi=100)
            plt.clf()
            return out_name

    def gen_simple_bar_chart(self,datax,datay,x_name,y_name,main_name,output_file_name,ret_obj=False):
        """generate a simple bar plot of a list of numbers"""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_title(main_name)
        ax.grid(True)
        out_name = output_file_name + ".png"
        plt.bar(datax,datay,color='blue')

        #set X axis limit
        plt.xlim(0,len(datax))
        
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:           
            plt.savefig(out_name,dpi=100)
            plt.clf()
            return out_name


    #for analyze_vcf.py.plot_variant_all_statuses_in_windows
    #generate a stacked bar chart. This requires a data structure that looks like
    #expects a list of np.array([])
    #where each item in the list is going to be a stack
    #number of items in the np.array defines how many columns there will be in the chart
    #note: all arrays have to be the same length
    #looks like:
    #[
    #np.array([1,2,3,4]),
    #no.array([4,5,6,7])
    #]
    #above will create 4 columns with 2 stacks each
    def gen_stacked_bar_chart(self,data,labels,y_label,main,out_name,vertical_separators_names,vertical_separators_pos,all_real_positions,ret_obj=False):
        """generate a stacked bar chart. This requires a data structure that looks like"""
        #print vertical_separators
        constant = PlatformConstant()
        num_cols=0
        fig = plt.figure(figsize=(24,12))
        ax = fig.add_subplot(111)
        try:
            num_cols = len(data[0])
        except:
            raise "ERROR: data set given to get_stacked_bar_chart empty"
        bar_locations = np.arange(num_cols)
        charts=[]
        running_total=data[0]
        #find global max of data set to print text at mid hieght
        global_max=0
        
        for i in range(0,len(data)):
            stats = SimpleStat(data[i])
            if stats.get_max()>global_max:
                global_max = stats.get_max()
            if i==0:
                charts.append(plt.bar(bar_locations, data[i],color=constant.PLOT_COLORS[i],edgecolor='none',label=labels[i]))
            else:
                charts.append(plt.bar(bar_locations, data[i],color=constant.PLOT_COLORS[i],bottom=running_total,edgecolor='none',label=labels[i]))
                running_total += data[i]
                
        if len(vertical_separators_names)>0:
            for i in range(0,len(vertical_separators_names)):
                name = vertical_separators_names[i]
                pos=vertical_separators_pos[i]
                ax.axvline(x=pos,color='k')
                if global_max/2 > 1:
                    ax.text(pos, -14, name  ,rotation=45)
                else:
                    ax.text(pos, 0.8, name  ,rotation=45)
                
        plt.ylabel(y_label)
        plt.title(main)
        plt.grid(False)

        #leg = plt.legend(charts,labels)
        #legend below axis
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),fancybox=True, shadow=True, ncol=5)

        from matplotlib.font_manager import FontProperties
        fontP = FontProperties()
        fontP.set_size('small')

        #set X axis limit
        plt.xlim(0,len(data[0]))

        #set ticks at every 40 places
        tick_locs=[]
        tick_labels=[]
        if len(data[0])>40:
            for i in range(0,len(all_real_positions)):
                pos = all_real_positions[i]
                if i %100 == 0:
                    tick_locs.append(i)
                    tick_labels.append(str(pos))
        else:
            for i in range(0,len(all_real_positions)):
                pos = all_real_positions[i]
                tick_locs.append(i)
                tick_labels.append(str(pos))


        plt.xticks(tick_locs, tick_labels)
        plt.savefig(out_name+'.png',dpi=100)



    #analyze_vcf.py.plot_snp_density_in_windows uses this plot
    def gen_multicolored_bar_chart(self,pass_counts_in_windows_by_chromosome,chrom_order,x_name,y_name,main_name,output_file_name,ret_obj=False):
        """generate a simple bar plot of a list of numbers"""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_title(main_name)
        ax.grid(True)
        out_name = output_file_name + ".png"

        prev_length=0
        total_length=0
        color_index=0
        xlabels=[]
        #find global max to figure out best place to put text labels
        global_max=0
        for chromosome_name,percentages in pass_counts_in_windows_by_chromosome.iteritems():
            stats = SimpleStat(percentages)
            if stats.get_max() > global_max:
                global_max = stats.get_max()

        #do the real plotting
        #for chromosome_name,percentages in pass_counts_in_windows_by_chromosome.iteritems():
        for chromosome_name in chrom_order:
            percentages = pass_counts_in_windows_by_chromosome[chromosome_name]
            xlabels.append(chromosome_name)
            if color_index % 2 == 0:
                ax.bar(range(prev_length,prev_length+len(percentages)),percentages, color="blue",edgecolor='none')
            else:
                ax.bar(range(prev_length,prev_length+len(percentages)),percentages, color="green",edgecolor='none')
            color_index+=1
            #draw labels now, using midpoints of global max and x position as coords for plotting text
            stats = SimpleStat(percentages)
            ax.text(prev_length+(len(percentages)/2), (global_max/2), chromosome_name  ,rotation=90)    
            total_length+=len(percentages)
            prev_length = prev_length+len(percentages)
            ax.axvline(x=prev_length,color='k')

        #set X axis limit
        plt.xlim(0,total_length)
        
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:           
            plt.savefig(out_name,dpi=100)
            plt.clf()
            return out_name



    #generates a bar plot, data=list of numbers. Please specify ouput
    #file name without any extensions.
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_bar_chart_with_axis(self,data,ls,x_name,y_name,main_name,output_file_name, ret_obj = False):
        """generate a simple bar plot of a list of numbers"""
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(main_name)
        labels = ax.set_xticklabels((ls))
        out_name = output_file_name + ".png"
        plt.ylabel(y_name)
        plt.xlabel(x_name)
        plt.grid(True)
        plt.plot(data)
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            plt.clf()
            return out_name



    #generates a histogram, data=list of numbers. Please specify ouput
    #file name without any extensions.
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_histogram(self,data,x_name,y_name,main_name,output_file_name, ret_obj=False):
        """generate a simple bar plot of a list of numbers"""
        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_title(main_name)
        ax.grid(True)
        out_name = output_file_name + ".png"
        
        #n, bins, patches = ax.hist(data, 50, normed=1, facecolor='blue', alpha=0.75)
        n, bins, patches = ax.hist(data, 50, facecolor='#0099FF', alpha=0.75)

        np.seterr(all='ignore')
        try:
            stats = SimpleStat(data)
            mu = stats.get_mean()
            sigma= stats.get_standard_deviation()
            sigma= stats.get_standard_deviation()
            bincs = 0.5*(bins[1:]+bins[:-1])
            y = mlab.normpdf( bincs, mu, sigma)
            l = ax.plot(bincs, y, 'r--', linewidth=1)
        except:
            print "plotting histogram: note: unable to calculate standard deviation and mean due to data characteristics"

        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:    
            plt.savefig(out_name,dpi=100)
            plt.clf()
            return out_name



    #1. analyze_vcf.py. plot_MQ_to_position uses this
    #2. identify_coverage_anomalies.py gen_scatter_plot function uses this
    #3. identify_coverage_anomalies.py  function coverage_delta_plots uses this (twice here)
    #generates a scatter plot, data_x=list of X, data_y=list of y.  Please specify ouput
    #file name without any extensions.
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_scatter_plot(self,
                         data_x,
                         data_y,
                         x_name,
                         y_name,
                         main_name,
                         output_file_name,
                         chrom_break_positions,
                         chrom_break_positions_labels,
                         all_real_positions=[],
                         ret_obj=False):
        """generate a simple scatter plot of a list of tuples"""
        if len(data_x) != len(data_y):
            raise "Error: lenghth of X " + str(len(data_x)) + " and length of Y " + str(len(data_y))
        out_name = output_file_name + ".png"
        matplotlib.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(24,6))
        ax = fig.add_subplot(111)
        ax.plot(data_x, data_y, '.',rasterized=True)

        #add chromosome breaks
        label_y_pos = max(data_y)/2
        label_y_pos = -12

        #chrom_break_positions,chrom_break_positions_labels,
        if len(chrom_break_positions)>0:
            for i in range(0,len(chrom_break_positions)):
                pos=chrom_break_positions[i]
                name = chrom_break_positions_labels[i]
                ax.axvline(x=pos,color='k')
                ax.text(pos,label_y_pos,name,fontsize=10,rotation=45)


        ax.set_title(main_name)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.grid=(True)

        #set X limit
        plt.xlim(0,max(data_x))

        #set ticks at every 40 places
        tick_locs=[]
        tick_labels=[]
        if len(all_real_positions)>40:
            for i in range(0,len(data_x)):
                pos = all_real_positions[i]
                if i %100 == 0:
                    tick_locs.append(i)
                    tick_labels.append(str(pos))
        if len(all_real_positions)<40 and len(all_real_positions) > 0:
            for i in range(0,len(data_x)):
                pos = all_real_positions[i]
                tick_locs.append(i)
                tick_labels.append(str(pos))
        plt.xticks(tick_locs, tick_labels)

        
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name



    #analyze_vcf.py.plot_variant_status uses this plot
    #plots muliple data sets on same plot
    #[[x1],[y1],[x2],[y2]]
    #here the pairs of indexs in input are considered x,y pair lists
    def gen_multiple_scatters_on_same_plot(self,all_data,x_name,y_name,main_name,legend,output_file_name,chrom_break_positions,chrom_break_positions_labels,ret_obj=False):
        """generate a multi scatter plot of a list of tuples"""
        #check if data ok
        for i in range(0,len(all_data),2):
            data_x=all_data[i]
            data_y=all_data[i+1]
            if len(data_x) != len(data_y):
                raise "Error: lenghth of X " + str(len(data_x)) + " and length of Y " + str(len(data_y))
        #plot now    
        out_name = output_file_name + ".png"
        matplotlib.rcParams['axes.unicode_minus'] = False
        fig = plt.figure(figsize=(24,12))
        ax = fig.add_subplot(111)

        colors=['r','g','y','b','p']
        c=0
        y_max=0
        for i in range(0,len(all_data),2):
            data_x=all_data[i]
            data_y=all_data[i+1]
            #find the global max to help with plotting labels
            stat = SimpleStat(data_y)
            if stat.get_max()>y_max:
                y_max=stat.get_max()
            plt.plot(data_x, data_y,'.',markersize=4,c=colors[c],rasterized=True, label=legend[c])
            x_max=max(data_x)
            c+=1

        #set x limit
        plt.xlim(0,x_max)

        #chrom_break_positions,
        #chrom_break_positions_labels,

        #add labels and chromosome breaks
        label_y_pos=y_max/2
        if len(chrom_break_positions)>0:
            #for name,pos in breaks.iteritems():
            for i in range(0,len(chrom_break_positions)):
                pos=chrom_break_positions[i]
                name = chrom_break_positions_labels[i]
                ax.axvline(x=pos,color='k')
                ax.text(pos,label_y_pos,name,fontsize=10,rotation='vertical')
            
        ax.set_title(main_name)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.grid=(True)

        from matplotlib.font_manager import FontProperties
        fontP = FontProperties()
        fontP.set_size('small')

        # Shink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])

        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),fancybox=True, shadow=True, ncol=10,prop=fontP)

        #plt.legend()
        #plt.legend(loc="upper left", bbox_to_anchor=(1,1))
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name


    
    #generates a mlticolored scatter plot,
    #data_x=list of X, data_y=list of y.  Please specify ouput
    def gen_multicolored_scatter_plot(self,data_x,data_y,names,x_name,y_name,main_name, output_file_name,ret_obj=False):
        """generate a simple scatter plot of a list of tuples"""
        if len(data_x) != len(data_y):
            raise "Error: lenghth of X " + str(len(data_x)) + " and length of Y " + str(len(data_y))
        out_name = output_file_name + ".png"
        matplotlib.rcParams['axes.unicode_minus'] = False
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(data_x, data_y, '.',rasterized=True)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_title(main_name)
        ax.grid=(True)

        #add names
        i=0
        for name in names:
            ax.text(data_x[i],data_y[i],name,fontsize=10)
            i+=1
            
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name    

    #takes in a data structure that looks like following
    #data={
    #      'column names':['col1','col2','col3','col4'],
    #      'ring1':[[20,22,24,28]],
    #      'ring2':[[18,14,32,24]],
    #      'ring3':[[32,34,28,22]],
    #      'ring4':[[12,12,20,20]]
    #      }
    #1. where all column names must be in the list under the key "column names"
    #2. There can be any number of rings
    #3. all values are "real" (IE not scaled via log or anything else)
    #last argument is a dictionary. 'grid':True turns on grid,
    #'col':<value> can feed in a specific color for color,
    #'ret_bj':True returns a plt object versus saving a png file
    #please contact harindra@broad for details if needed.
    def gen_radar_plot(self,data,axis_count,main_name,out_name,options):
        """generate a radar plot"""
        radar_helper= RadarChart(axis_count)
        theta = radar_helper.radar_factory()
        spoke_labels = data.pop('column names')

        fig = plt.figure(figsize=(8 ,8))
        fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)
        ax = fig.add_subplot(1, 1, 1, projection='radar')
        #add title
        ax.set_title(main_name, weight='bold', size='medium', position=(0.5, 1.1),horizontalalignment='center', verticalalignment='center')
        colors=['r','b','g','y','c','m','k']
        color_index=0
        for n, title in enumerate(data.keys()):
            for d in data[title]:
                if options.has_key('col') and options['col'] == '-1':
                    ax.plot(theta, d, color=colors[color_index],rasterized=True,label=title)
                else:
                    ax.plot(theta, d, color=col,rasterized=True,label=title)
                color_index+=1
        ax.set_varlabels(spoke_labels)
        if options.has_key('grid'):
            ax.grid(options['grid'])
        else:
            ax.grid(False)
        plt.legend()
        
        if options.has_key('ret_obj') and options['ret_obj']:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name+".png",facecolor='w', edgecolor='none')
            return out_name


    def merge_2_image_charts(self,top,bottom,out_name):
        """take in two images, make top main chart and bottom the watermark"""
        im_top = image.imread(top)
        im_bottom = image.imread(bottom)
        im_bottom[:,:,-1] = 0.4
        fig = plt.figure(1,(10,10))
        ax = fig.add_subplot(111)
        ax.set_visible(False)
        fig.figimage(im_top,0.5,0.5)
        fig.figimage(im_bottom,0.5,0.5)
        plt.savefig(out_name+".basecomp_variant_eyeplots.merged.png", edgecolor='none')


    #generate a heat map
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_heat_map(self,data,x_name,main_title,colorbar_label,output_file_name,ret_obj=False):
        """to generate a heat map"""
        out_name = output_file_name + ".png"        
        fig = plt.figure()
        main = plt.imshow(data, aspect='auto',cmap='Blues')
        main.set_label(main_title)
        cb = plt.colorbar()
        cb.set_label(colorbar_label)

        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name




    #generate a heat map
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_heat_map_advanced(self,data,x_name,main_title,colorbar_label,output_file_name,ret_obj=False):
        """to generate a heat map"""
        out_name = output_file_name + ".png"        
        fig = plt.figure()
        main = plt.imshow(data, aspect='auto',cmap='Reds',interpolation='nearest')
        main.set_label(main_title)
        cb = plt.colorbar()
        cb.set_label(colorbar_label)

        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name

    #following tools use this:
    # 1. analyze_vcf_file.py plot_distribution_of_pass_fail_states function
    #
    #generate pie plot. Data = list of numbers, list of labels
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_pie_plot(self,data,lbls,main_title,output_file_name,ret_obj=False):
        """to generate a pie plot"""
        out_name = output_file_name + ".png"
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111)
        ax.set_title(main_title)
        patches, texts, autotexts = plt.pie(data, explode=None, labels=lbls,autopct='%1.1f%%', pctdistance=0.6, labeldistance=1.025, shadow=False)

        proptease = fm.FontProperties()
        proptease.set_size('small')
        plt.setp(autotexts, fontproperties=proptease)
        plt.setp(texts, fontproperties=proptease)
        if ret_obj:
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name


    #composite chart function. Takes in a list of list of numbers. list[list]. Here each of inner lists
    #represent a single plot
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_composite_plot(self,all_plot_values,x_name,y_name,main_name,output_file_name,ret_obj=False):
       """plots a list of, list of numbers in the same chart for example: list[list]
          each i is a separate plot, and each j of the i is a value for a plot """
       fig = plt.figure()
       ax = fig.add_subplot(111)
       ax.set_xlabel(x_name)
       ax.set_ylabel(y_name)
       ax.set_title(main_name)
       ax.grid(True)
       for single_plot_values in all_plot_values:
           plt.plot(single_plot_values)
           
       
       out_name = output_file_name + ".png"
       if ret_obj:
           chrt = Chart(plt,fig)
           return chrt
       else:
           plt.save_fig(out_name,dpi=100)
           return out_name


    #ACTIVE
    #composite chart function. Takes in a list of list of numbers. list[list]. Here each of inner lists
    #represent a single plot
    #if a "True" is given for ret_obj arg, a Chart object is returned versus a file
    #file being saved according to the name given. By default, a file is saved. Use the
    #file_name to specify where the file is to be saved. IE give full path to file name.
    def gen_line_plot_with_two_y_scales(self,smaller_plot_values,
                                        smaller_plot_labels,
                                        x_name,
                                        y1_name,
                                        large_plot_values,
                                        large_plot_label,
                                        y2_name,
                                        main_name,
                                        output_file_name,
                                        scaf_break_coords,
                                        scaf_break_labels,
                                        ret_obj=False,
                                        both_xy_coords_given=False):
       """plots a list of, list of numbers in the same chart for example: list[list]
          each i is a separate plot, and each j of the i is a value for a plot """
       fig = plt.figure(figsize=(24,6))
       ax = fig.add_subplot(111)
       ax.set_xlabel(x_name)
       ax.set_ylabel(y1_name)
       ax.set_title(main_name)
       ax.grid(True)
       index=0
       single_plot_values=[]
       colors=('green','blue','yellow','magenta','purple','gray','black')
       while index < len(smaller_plot_values):
           single_plot_values = smaller_plot_values[index]
           if both_xy_coords_given==False:
                plt.plot(single_plot_values,label=smaller_plot_labels[index],color=colors[index])
           else:
                plt.plot(single_plot_values[0],single_plot_values[1],label=smaller_plot_labels[index],color=colors[index])
           index+=1

       #add scaffold breaks
       scaf_index=0
       for coord in scaf_break_coords:
           ax.axvline(x=coord,color='k')
           ax.text(coord,-3,scaf_break_labels[scaf_index],fontsize=10,rotation='vertical')
           scaf_index+=1

       from matplotlib.font_manager import FontProperties
       fontP = FontProperties()
       fontP.set_size('small')
       
       #Srhink current axis's height by 10% on the bottom
       box = ax.get_position()
       ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])

       # Put a legend below current axis
       ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),fancybox=True, shadow=True, ncol=10,prop=fontP)
       
       other_y_axis = ax.twinx()
       if both_xy_coords_given==False:
           other_y_axis.plot(large_plot_values, label=large_plot_label, color='red')
       else:
           other_y_axis.plot(large_plot_values[0],large_plot_values[1], label=large_plot_label, color='red')

       other_y_axis.set_ylabel(y2_name)
       out_name = output_file_name + ".png"

       #set x limit
       if both_xy_coords_given==True:
           plt.xlim(0,len(single_plot_values[0]))
       else:
           plt.xlim(0,len(single_plot_values))
       
       if ret_obj:
           chrt = Chart(plt,fig)
           return chrt
       else:
           current_fig = plt.gcf()
           current_fig.savefig(out_name,dpi=100)
           return out_name


    #rotates a list that looks like  [ [1],[2],[3],[4] ] to look like [ [1,2,3,4] ]
    def rotate_2d_array(self,twod):
        len_of_orig_row = len(twod[0])
        rotated=[]
        for i in range(0,len_of_orig_row):
            rotated.append([])
        for row in twod:
            for col in range(0,len(row)):
                rotated_row = col
                rotated[rotated_row].append(row[col])
        return rotated








   
    #ACTIVE
    def gen_line_plot_with_external_gc_plot(self,smaller_plot_values,smaller_plot_labels,x_name,y1_name,large_plot_values,large_plot_label,y2_name,main_name,output_file_name,scaf_break_coords,scaf_break_labels,ret_obj=False):
       """plots a list of, list of numbers in the same chart for example: list[list]
          each i is a separate plot, and each j of the i is a value for a plot """
       fig = plt.figure(figsize=(24,18))
       fontP = FontProperties()
       fontP.set_size('small')

       #add GC% plot as a subplot
       other_y_axis = fig.add_subplot(211)
       other_y_axis.set_title(main_name)
       other_y_axis.plot(large_plot_values, label=large_plot_label, color='red')

       #make gc% plot x axis and ticks invisible
       #other_y_axis.set_frame_on(False)
       #plt.setp(other_y_axis.get_xticklabels(), visible=False)
       
       other_y_axis.set_ylabel(y2_name)
       out_name = output_file_name + ".png"

       #set x limit
       plt.xlim(0,len(smaller_plot_values[0]))

       #add scaffold breaks
       scaf_index=0
       for coord in scaf_break_coords:
           other_y_axis.axvline(x=coord,color='k')
           other_y_axis.text(coord,-3,scaf_break_labels[scaf_index],fontsize=10,rotation=45)
           scaf_index+=1
       
       ax = fig.add_subplot(212)
       ax.set_xlabel(x_name)
       ax.set_ylabel(y1_name)
       #ax.grid(True)
       
       # Put a legend below current axis
       #other_y_axis.legend(loc=(10,10),bbox_to_anchor=(0.02, -0.15),fancybox=True, shadow=True, ncol=10,prop=fontP)

       
       index=0
       single_plot_values=[]
       colors=('green','blue','yellow','magenta','purple','gray','black')
       while index < len(smaller_plot_values):
           single_plot_values = smaller_plot_values[index]
           plt.plot(single_plot_values,label=smaller_plot_labels[index],color=colors[index])
           index+=1
           
       #set x limit
       plt.xlim(0,len(single_plot_values))

       #add scaffold breaks
       scaf_index=0
       for coord in scaf_break_coords:
           ax.axvline(x=coord,color='k')
           ax.text(coord,-3.4,scaf_break_labels[scaf_index],fontsize=10,rotation=45)
           scaf_index+=1
       
       #Srhink current axis's height by 10% on the bottom
       box = ax.get_position()
       ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])

       # Put a legend below current axis
       ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),fancybox=True, shadow=True, ncol=10,prop=fontP)
       
       
       if ret_obj:
           chrt = Chart(plt,fig)
           return chrt
       else:
           plt.save_fig(out_name,dpi=100)
           return out_name



#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


   
    #ACTIVE
    def gen_line_plot_with_layout_managed_gc_plot(self,smaller_plot_values,smaller_plot_labels,x_name,y1_name,large_plot_values,large_plot_label,y2_name,main_name,output_file_name,scaf_break_coords,scaf_break_labels,minor_xtick_coords, minor_xtick_labels,ret_obj=False):
       """plots a list of, list of numbers in the same chart for example: list[list]
          each i is a separate plot, and each j of the i is a value for a plot """
       fig = plt.figure(figsize=(24,12))
       subplots = gridspec.GridSpec(2, 1, height_ratios=[1, 3]) 
       fontP = FontProperties()
       fontP.set_size('small')

       #add GC% plot as a subplot
       other_y_axis = fig.add_subplot(subplots[0])
       
       other_y_axis.set_title(main_name)
       other_y_axis.plot(large_plot_values, label=large_plot_label, color='red')

       #ticks
       other_y_axis.set_xticks([])

       #make gc% plot x axis and ticks invisible
       plt.setp(other_y_axis.get_xticklabels(), visible=False)

       #gc plot axis labels
       other_y_axis.set_ylabel(y2_name)

       out_name = output_file_name + ".png"

       #set x limit
       plt.xlim(0,len(smaller_plot_values[0]))

       #add scaffold breaks
       scaf_index=0
       for coord in scaf_break_coords:
           other_y_axis.axvline(x=coord,color='k')
           scaf_index+=1

       #add y limit
       plt.ylim(0,100)
       

       #add legend
       plt.legend(loc="upper left")

       #add next plot
       ax = fig.add_subplot(subplots[1])
       ax.set_xlabel(x_name,labelpad=60)
       ax.set_ylabel(y1_name)

       #only add ticks to left side of chart, take out right, and top
       ax.yaxis.tick_left()
       ax.xaxis.tick_bottom()
              
       index=0
       single_plot_values=[]
       colors=('green','blue','yellow','magenta','purple','gray','black')
       global_y_max=0
       global_avrg_max = 0
       global_stddev_max=0
       while index < len(smaller_plot_values):
           single_plot_values = smaller_plot_values[index]
           plt.plot(single_plot_values,label=smaller_plot_labels[index],color=colors[index])

           #gather some stats to check if ylim is needed
           stats = SimpleStat(single_plot_values)
           local_max = stats.get_max()
           local_avrg_max = stats.get_mean()
           local_stddev_max = stats.get_standard_deviation()
           if local_max>global_y_max:
               global_y_max = local_max
           if local_avrg_max > global_avrg_max:
               global_avrg_max = local_avrg_max
           if local_stddev_max > global_stddev_max:
               global_stddev_max = local_stddev_max

           index+=1
           
       #set x limit
       plt.xlim(0,len(single_plot_values))

       #see if ylim is needed (if max is more than 4xstdevs)
       if global_y_max > (global_avrg_max + (4 * global_stddev_max)):
           plt.ylim(0,(global_avrg_max + (4 * global_stddev_max)))

       #add scaffold breaks
       scaf_index=0
       for coord in scaf_break_coords:
           ax.axvline(x=coord,color='k')
           if global_y_max>100:
                ax.text(coord,-18.0,scaf_break_labels[scaf_index],fontsize=10,rotation=90)
           else:
                ax.text(coord,-2.0,scaf_break_labels[scaf_index],fontsize=10,rotation=90)
           scaf_index+=1

       #add minor xticks and labels
       ax.xaxis.set_ticks(minor_xtick_coords)
       for i in range(0,len(minor_xtick_coords)):
           ax.text(minor_xtick_coords[i],-12.0,minor_xtick_labels[i],fontsize=8,rotation=90)
       
       #bring plots closer together
       plt.subplots_adjust(hspace = .042)

       #take out default x labels
       ax.set_xticklabels(())

       #legend below axis
       box = ax.get_position()
       ax.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
       ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),fancybox=True, shadow=True, ncol=5)
       
       if ret_obj:
           chrt = Chart(plt,fig)
           return chrt
       else:
           plt.save_fig(out_name,dpi=100)
           return out_name




#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,



    #data is a [], list[0]=GC,  list[1:]=each is a dict. with values being a list of cvrg from a bam. key
    #is bam name.
    #xlim is real value that should be the xlimit of x axis
    #ylim_multiple is the number std should be multiplied to get the ylim of y axis (this
    #is to reduce effects of outliers
    def gen_line_and_heatmap_in_same_chart(self,all_chart_data,main,xlim=-1,ylim_multiple=-1):
        fig = plt.figure(figsize=(12, 12))
        fig.subplots_adjust(hspace=2)
        
        num_cols=1
        num_rows = len(all_chart_data[1])+1
        index=0
        axes=[]

        #heat map
        single=all_chart_data[0]
        single_chart_data = self.rotate_2d_array(single)
        axes.append(fig.add_subplot(num_rows,num_cols,0))
        axes[index].imshow(single_chart_data, aspect='auto',interpolation='nearest',cmap='Blues')
        axes[index].set_title("GC% of Sequence Windows")
        axes[index].set_xlabel("sequence window")
        
        #line charts
        num_fasta_blocks=0
        cvrg_charts = all_chart_data[1]
        index=1
        single_chart_data=[]
        x_coordinate=0
        scaffold_break_coords=[]
        coords=[]
        for bam_name,values in cvrg_charts.iteritems():
            scaffold_break_coords=[]
            scaffold_labels=[]
            for fasta_id,cvrg_values in values.iteritems():
                num_fasta_blocks+=1
                for single_cvrg in cvrg_values:
                    single_chart_data.append(single_cvrg)
                    coords.append(x_coordinate)
                    x_coordinate+=1
                scaffold_break_coords.append(x_coordinate)
                scaffold_labels.append(fasta_id)


            axes.append(fig.add_subplot(num_rows,num_cols,index))            
            #line plot
            axes[index].plot(single_chart_data)
            axes[index].set_title(bam_name + main + " line plot")
            axes[index].set_xlabel("sequence window")
            axes[index].set_ylabel("coverage")
            plt.minorticks_on()
            
            #testing of adding further granularity to minor ticks
            #from matplotlib.ticker import MultipleLocator, AutoLocator, FixedFormatter
            #axes[index].xaxis.set_major_locator(AutoLocator()) 
            #major = axes[index].xaxis.get_majorticklocs() 
            #minor =  (major[-1]-major[0])/(len(major)-1) /10.
            #minorticklabels = []
            #i=0
            #for x in coords:
            #    if i%400==0:
            #        minorticklabels.append(x)
            #    i+=1
            #axes[index].xaxis.set_minor_formatter(FixedFormatter(minorticklabels))
            #axes[index].xaxis.set_minor_locator(MultipleLocator(minor))

            #add scaffold breaks
            for scaffold_break_coord in scaffold_break_coords:
                axes[index].axvline(x=scaffold_break_coord,color='k')
            #label X-axis (if enough fasta blocks found)
            if num_fasta_blocks>1:
                locs,labels = plt.xticks(scaffold_break_coords,scaffold_labels,rotation='vertical')
            else:
                print "could not activate ticks (probably only single fasta)skipping this chart feature,number of fasta blocks found: "+str(num_fasta_blocks)
            
            index+=1


        #set limits
        if xlim != -1:
            axes[index].xlim(0,xlim)
        if xlim == -1:
            plt.xlim(0,len(single_chart_data))
                    
        if ylim_multiple != -1:
            stats = SimpleStat(single_chart_data)
            ylim= (stats.get_standard_deviation() * ylim_multiple)+stats.get_mean()
            plt.ylim(0,ylim)            
        chrt = Chart(plt,fig)
        return chrt




    #all_cgart_data[0] is the data for the GC% heat map. This ALWAYS should be the heat map in the format,
    #[[]]
    def gen_line_and_heatmap_in_same_chart_integrated(self,all_chart_data,main,xlim=-1,ylim_multiple=-1):
        fig = plt.figure(figsize=(24, 6))
        fig.subplots_adjust(hspace=2)                
        num_cols=1
        num_rows = len(all_chart_data[1])
        index=0
        axes=[]

        #print all_chart_data
   
        #first data items is GC% and meant for heat map. First change it to s 2x1 array for heat map
        single=all_chart_data[0]
        rotated_chart_data = self.rotate_2d_array(single)

        #next data item is a collection of coverage data, one per BAM
        cvrg_charts = all_chart_data[1]
        
        #go through each bam coverage collection and generate BAM plots with GC% in it
        single_chart_data=[]
        for bam_name,values in cvrg_charts.iteritems():
            #some counters to track scaffold breaks and their coordinates
            num_fasta_blocks=0
            x_coordinate=0
            coords=[]

            #fasta specific vars
            scaffold_break_coords=[]
            scaffold_labels=[]

            #sort keys (scaffolds are named in ascending order (scaffold01..02 etc)
            value_keys = values.keys()
            value_keys.sort()
            for fasta_id in value_keys:
                cvrg_values = values.get(fasta_id)
                print fasta_id,str(len(cvrg_values))
                num_fasta_blocks+=1
                for single_cvrg in cvrg_values:
                    single_chart_data.append(single_cvrg)
                    coords.append(x_coordinate)
                    x_coordinate+=1
                scaffold_break_coords.append(x_coordinate)
                scaffold_labels.append(fasta_id)
                    
            #line plot
            axes.append(fig.add_subplot(num_rows,num_cols,index))
            axes[index].plot(single_chart_data)
            axes[index].set_title(bam_name + main + " line plot")
            axes[index].set_xlabel("sequence window")
            axes[index].set_ylabel("coverage")
            plt.minorticks_on()

            #set limits
            if xlim != -1:
                axes[index].xlim(0,xlim)
            if xlim == -1:
                plt.xlim(0,len(single_chart_data))

            if ylim_multiple != -1:
                stats = SimpleStat(single_chart_data)
                ylim= (stats.get_standard_deviation() * ylim_multiple)+stats.get_mean()
                plt.ylim(0,ylim)


            #add scaffold breaks
            for scaffold_break_coord in scaffold_break_coords:
                axes[index].axvline(x=scaffold_break_coord,color='k')
            #label X-axis (if enough fasta blocks found)
            if num_fasta_blocks>1:
                locs,labels = plt.xticks(scaffold_break_coords,scaffold_labels,rotation='vertical')
            else:
                print "could not activate ticks (probably only single fasta)skipping this chart feature,number of fasta blocks found: "+str(num_fasta_blocks)
            

            #add heat map
            axes[index].set_axes(plt.axes([0.126,0.79,0.772,0.1]))
            axes[index].get_axes().imshow(rotated_chart_data, aspect='auto',interpolation='nearest',cmap='Blues')
            axes[index].get_axes().get_yaxis().set_visible(False)
            axes[index].get_axes().get_xaxis().set_visible(False)


            index+=1
            
        chrt = Chart(plt,fig)
        return chrt



    #ACTIVE
    #take #3
    #-------------------
    #data_for_gc_heatmap  : [[gc1], [gc2],[gc3], [gc4] .... ]
    #data_for_line_plot   : 'bam_name': { 'fasta_id1':[cvr1,cvr2...], 'fasta_id2':[cvrg1,cvrg2]  }
    def gen_line_and_heatmap_in_same_chart_integrated_single(self,data_for_gc_heatmap, data_for_line_plot,main,xlim=-1,ylim_multiple=-1):
        fig = plt.figure(figsize=(24, 12))
        fig.subplots_adjust(hspace=2)
        axe = fig.add_subplot(111)

        #First change heatmap data to a 2x1 array for heat map
        rotated_chart_data = self.rotate_2d_array(data_for_gc_heatmap)

        #go through each bam coverage collection and generate BAM plots with GC% in it
        single_chart_data=[]
        for bam_name,values in data_for_line_plot.iteritems():
            #some counters to track scaffold breaks and their coordinates
            num_fasta_blocks=0
            x_coordinate=0
            coords=[]

            #fasta specific vars
            scaffold_break_coords=[]
            scaffold_labels=[]

            #sort keys (scaffolds are named in ascending order (scaffold01..02 etc)
            value_keys = values.keys()
            value_keys.sort()
            for fasta_id in value_keys:
                cvrg_values = values.get(fasta_id)
                num_fasta_blocks+=1
                for single_cvrg in cvrg_values:
                    single_chart_data.append(single_cvrg)
                    coords.append(x_coordinate)
                    x_coordinate+=1
                scaffold_break_coords.append(x_coordinate)
                scaffold_labels.append(fasta_id)

            #line plot
            axe.plot(single_chart_data)
            axe.set_title(bam_name + main + " line plot")
            axe.set_xlabel("sequence window")
            axe.set_ylabel("coverage")
            plt.minorticks_on()

            #set limits
            if xlim != -1:
                axes.xlim(0,xlim)
            if xlim == -1:
                plt.xlim(0,len(single_chart_data))

            if ylim_multiple != -1:
                stats = SimpleStat(single_chart_data)
                ylim= (stats.get_standard_deviation() * ylim_multiple)+stats.get_mean()
                plt.ylim(0,ylim)


            #add scaffold breaks
            for scaffold_break_coord in scaffold_break_coords:
                axe.axvline(x=scaffold_break_coord,color='k')
            #label X-axis (if enough fasta blocks found)
            if num_fasta_blocks>1:
                locs,labels = plt.xticks(scaffold_break_coords,scaffold_labels,rotation='vertical')
            else:
                print "could not activate ticks (probably only single fasta)skipping this chart feature,number of fasta blocks found: "+str(num_fasta_blocks)


            #add heat map
            axe.set_axes(plt.axes([0.126,0.79,0.965,0.112]))
            heatmap_axe= axe.get_axes().imshow(rotated_chart_data, aspect='auto',interpolation='nearest',cmap='Blues')
            axe.get_axes().get_yaxis().set_visible(False)
            axe.get_axes().get_xaxis().set_visible(False)
            plt.colorbar(heatmap_axe, ticks=np.linspace(0,100))

        chrt = Chart(plt,fig)
        return chrt



    #primary user:read analysis tool
    #all_chart_data = { 'legend':([],[]), 'legend':([],[]) }
    def gen_multiple_line_plots(self,all_chart_data,main,x_label,y_label,out_name,options=dict()):
        fig = plt.figure()
        fig.subplots_adjust(hspace=2)

        axe = fig.add_subplot(111)  
        for name,values in all_chart_data.iteritems():
            axe.plot(values[0],values[1], label=name)
            
        single_chart_data=values
        #set limits
        if options.has_key('xlim') and  options['xlim'] != 0:
            axe.xlim(0,xlim)
            
        if options.has_key('xlim') and  options['xlim'] == 0:
            plt.xlim(0,len(single_chart_data))
                    
        if options.has_key('ylim_multiple'):
            ylim_multiple = options['ylim_multiple']
            stats = SimpleStat(single_chart_data)
            ylim= (stats.get_standard_deviation() * ylim_multiple)+stats.get_mean()
            plt.ylim(0,ylim)

        axe.set_title(main)
        axe.set_xlabel(x_label)
        axe.set_ylabel(y_label)
        #plt.legend()

        from matplotlib.font_manager import FontProperties
        fontP = FontProperties()
        fontP.set_size('xx-small')

        # Shink current axis's height by 10% on the bottom
        box = axe.get_position()
        axe.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])

        # Put a legend below current axis
        axe.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),fancybox=True, shadow=True, ncol=1,prop=fontP)

        if options.has_key('ret_obj'):
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name

    

    #all_chart_data = { 'legend':([],[]), 'legend':([],[]) }
    def gen_multiple_scatter_plots(self,all_chart_data,main,x_label,y_label,out_name,options=dict()):
        fig = plt.figure()
        fig.subplots_adjust(hspace=2)
        
        axe = fig.add_subplot(111)  
        colors=['red','blue','green','cyan','magenta','black']
        markers=['+','x','d','^','>','v']
        index=0
        for name,values in all_chart_data.iteritems():
            axe.scatter(values[0],values[1], label=name,color=colors[index], s=4,marker=markers[index])
            index+=1
        
        axe.grid(True)
        plt.ylim(0,40)
        
        single_chart_data=values
        #set limits
        if options.has_key('xlim') and  options['xlim'] != 0:
            axe.xlim(0,xlim)
        
        if options.has_key('xlim') and  options['xlim'] == 0:
            plt.xlim(0,len(single_chart_data))
        
        if options.has_key('ylim_multiple'):
            ylim_multiple = options['ylim_multiple']
            stats = SimpleStat(single_chart_data)
            ylim= (stats.get_standard_deviation() * ylim_multiple)+stats.get_mean()
            plt.ylim(0,ylim)
        
        axe.set_title(main)
        axe.set_xlabel(x_label)
        axe.set_ylabel(y_label)
        #plt.legend()
        
        from matplotlib.font_manager import FontProperties
        fontP = FontProperties()
        fontP.set_size('small')
        
        # Shink current axis's height by 10% on the bottom
        box = axe.get_position()
        axe.set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
        
        # Put a legend below current axis
        axe.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),fancybox=True, shadow=True, ncol=1,prop=fontP)
        
        if options.has_key('ret_obj'):
            chrt = Chart(plt,fig)
            return chrt
        else:
            plt.savefig(out_name,dpi=100)
            return out_name

    

    #generate multiple histograms in same chart file
    def gen_histograms_in_same_chart(self,coverage_by_bam):
        fig = plt.figure()
        fig.subplots_adjust(hspace=2)
        
        num_cols=1
        num_rows = len(coverage_by_bam)

        index=0
        axes=[]
        for bam_name,fasta_cvrgs in coverage_by_bam.iteritems():
            all_cvrgs_from_single_bam=[]
            for fasta_name,cvrgs in fasta_cvrgs.iteritems():
                for cvrg in cvrgs:
                    all_cvrgs_from_single_bam.append(float(cvrg))
            #histogram(s)
            axes.append(fig.add_subplot(num_rows,num_cols,index))
            n, bins, patches = axes[index].hist(all_cvrgs_from_single_bam, 200, facecolor='blue', alpha=0.75)
            bincs = 0.5*(bins[1:]+bins[:-1])

            stats = SimpleStat(all_cvrgs_from_single_bam)
            mu = stats.get_mean()
            sigma= stats.get_standard_deviation()
            
            y = mlab.normpdf( bincs, mu, sigma)
            axes[index].plot(bincs, y, 'r--', linewidth=1)
            axes[index].set_title(bam_name + " sequence window histogram")
            axes[index].set_xlabel(" coverage bin ")
            axes[index].set_ylabel(" frequency ")            
            index+=1
            
        chrt = Chart(plt,fig)
        return chrt

    def surface_plot(self,x,y,z,out_name):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')        
        ax.plot_wireframe(x, y, z, rstride=10, cstride=10)
        plt.savefig(out_name,dpi=100)

#----tests
#cu = ChartUtilities()



#data_for_gc_heatmap  : 'fasta_name : [gc1, gc2, .... ]
#data_for_line_plot   : 'bam_name': { 'fasta_id1':[cvr1,cvr2...], 'fasta_id2':[cvrg1,cvrg2]  }

#line_data = {'bamname1':{'fasta_id1':[12,10,12,10,8,14,8,11,12,20,18],
#                         'fasta_id2':[10,4,6,8,8,11,18,17,11,19,18]
#                        }
#            }
#gc_data = [[40],[40],[39],[38],[40],[44],[40],[48],[38],[40],[44],[40],[40],[39],[38],[40],[44],[40],[48],[38],[40],[44]]
#chrt = cu.gen_line_and_heatmap_in_same_chart_integrated_single(gc_data, line_data,'main')
#chrt.save_as('line_and_gc_together.test.png')

#multi plot
#data1 = [1,2,3,4]
#data2 = [4,5,6,7]
#data3 = [6,7,8,9]
#data4 = [10,11,12,13,14,14,16]
#all_data = [data1,data2,data3,data4,data5]
#obj = cu.gen_composite_plot(all_data, "xname" , "yname" , "mainname" , "output_all_name" , True)
#obj.save_as("composite_test.png")

#composite plot2 (line with 2 scales/ys)
#data1 = [1,2,3,4]
#data2 = [4,5,6,7]
#data3 = [6,7,8,9]
#data4 = [10,11,12,13,14,14,16]
#large_data_set = [400,200,250,280,422,444,324]
#smaller_plot_values = [data1,data2,data3,data4]
#obj = cu.gen_line_plot_with_two_y_scales(smaller_plot_values,['data1','data2','data3','data4'],"x_name","y1_name",large_data_set,'large',"y2_name","main_name","output_file_name",True)
#obj.save_as("multi_axis_test.png")


#external GC% with line plotd
#data1 = [1,2,3,4,5,6,7]
#data2 = [4,5,6,7,3,2,5]
#data3 = [6,7,8,9,2,8,7]
#data4 = [10,11,12,13,14,14,16]
#gc = [40,20,25,28,42,44,32]
#lines = [data1,data2,data3,data4]
#obj = cu.gen_line_plot_with_external_gc_plot(lines,['data1','data2','data3','data4'],'x_name','y1_name',gc,'gc_label','y2_name','main_name','multiline_w_ext_gc',[4],['s1'],True)

#obj = cu.gen_line_plot_with_layout_managed_gc_plot(lines,['data1','data2','data3','data4'],'x_name','y1_name',gc,'gc_label','y2_name','main_name','multiline_w_ext_gc',[4],['s1'],[1,4,6,5],['1','4','6','5'],True)

#obj.save_as("ext_gc_with_lines.png")







#bar chart
#data = [1,2,3,4]
#obj = cu.gen_bar_chart(data,"xname","yname","mainname","output_file_name",True)
#print obj.__str__()


#scatter
#x = [1,2,3,4]
#y = [1,2,3,5]
#cu.gen_scatter_plot(x,y,"x","y","m","out")

#heat
#cu.gen_heat_map(np.clip(randn(250, 250), -1, 1),"x","main","out")

#pie
#cu.gen_pie_plot([10,20,30,40],["10","20","30","40"],"Pie","out")


#x=[[1,2,3]]
#y=[[4,6,8]]
#z=[[10,12,18]]
#cu.surface_plot(x,y,z"out")


#hist test
#ue_cn = [12,14,12,12,123,123,12,112,12,12,14,0,10,122,122,110,111,121,124]
#cu.gen_histogram(ue_cn,"Test X name","Test Y name","Test main name","Test_hist")

#radar test
#data={
#    'column names':['col1','col2','col3','col4'],
#    'ring1':[[20,22,24,28]],
#    'ring2':[[18,14,32,24]],
#    'ring3':[[32,34,28,22]],
#    'ring4':[[12,12,20,20]]
#    }
#xtraoptions={'grid':False, 'ret_obj':False, 'col':'-1'}
#rd = cu.gen_radar_plot(data,4,"main name","out_name_radar",xtraoptions)


#stacked chart test
#data = [np.array([1,2,3,4]),
#        np.array([10,12,14,18]),
#        np.array([22,14,11,11])
#        ]
#cu.gen_stacked_bar_chart(data,"y_label","main","stacked_test")



#test gen_multiple_line_plots
#c1=([1,2,3,4],[111,200,80,111])
#c2=([4,22,1,11],[121,124,111,112])
#c3=([1,2,3,4],[121,124,112,114])
#c4=([2,4,6,8],[119,114,116,118])
#all_chart_data= {'c1':c1, 'c2':c2, 'c3':c3, 'c4':c4}
#cu.gen_multiple_line_plots(all_chart_data,"main","x_label","y_label","test_out_name")
