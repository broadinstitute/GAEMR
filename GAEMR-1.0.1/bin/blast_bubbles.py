#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

from optparse import OptionParser
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot, mpl
from pylab import *
from scipy import *
from scipy import stats
from matplotlib.font_manager import FontProperties
from operator import itemgetter, attrgetter
import re

parser = OptionParser(usage="usage: %prog [options] <detail.contig.table.txt> <taxonomy.txt>")

parser.add_option('--output', '-o', action="store", default="blast_bubbles", type='string', dest="output",
                  help='Parsed blast output file header (default=%default).  Will write output.png.')

parser.add_option('--taxlevel', '-t', action="store", default="genus", type='string', dest="taxlevel",
                  help='Taxonomic level at which to color bubles (default=%default)')

parser.add_option('--delimiter', '-d', action="store", default=constant.TABLE_DELIMITER, dest="delimiter",
                  help='Contig detail table delimiter (default=%default)')

parser.add_option('--plotall', '-p', action="store_true", default=None, dest="plotall",
                  help='Option to plot all coverages.  No autoscaling.')

parser.add_option('--details', '-v', action="store", default=None, dest="details",
                  help='Output file header for plotting details (default=%default)')
(options, args) = parser.parse_args()

if len(args) < 2:
    parser.error("Must supply detailed contig table and taxonomy output.  These two files must be produced from the same fasta sequence.")
    sys.exit(-1)


def main():

    details_file = ''
    details_output = sys.stdout
    if options.details:
        details_file = options.details + '.blast_bubble_plot_details.out'
        try:
            details_output = open(details_file,'w')
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1

    details_output.write("Opening contig info details file: {}\n".format(args[0]))
    unitdata = {}
    for line in open(args[0], 'r'):
        if line.startswith("#"):
            continue
        unit, group, length,gc,cov_str,org,id = line.rstrip().split(options.delimiter)
        unitdata[unit] = []
        unitdata[unit].append(int(re.sub(",","",length)))
        unitdata[unit].append(gc)
        cov,rest = cov_str.split(" ")
        unitdata[unit].append(int(re.sub(",","",cov)))
        unitdata[unit].append(id)

    organisms = set()
    for line in open(args[1], 'r'):

        if line.startswith("#"):
            continue
    
        ctg,ctglen,hitlen,percov,tax = line.rstrip().split(options.delimiter)

        if (ctg != "#QueryId"):

            thisdata = {}

            for tax_info in tax.split(";"):
                #Hack to overcome issue with = in strain name
                hack = []
                hack = tax_info.split("=")
                #cls,val = tax_info.split("=")
                cls = hack[0]
                val = hack[1]
                thisdata[cls] = val
                
            if (options.taxlevel in thisdata.keys()):
                unitdata[ctg].append(thisdata[options.taxlevel])
                organisms.add(thisdata[options.taxlevel])
            else:
                unitdata[ctg].append(thisdata["domain"])
                organisms.add(thisdata["domain"])
        

    organisms = list(organisms)
    gcs = []
    coverages = []
    lengths = []
    colors = []
    nhgcs = []
    nhcoverages = []
    nhlengths = []
    org_count = {}


    cmap = cm.gist_rainbow
    num_colors = len(organisms)+1
    
    sm = cm.ScalarMappable()
    sm.set_cmap(cmap)
    sm.set_array(arange(0,num_colors-1,1))
    sm.autoscale()
    
    for ctg,data in sorted(unitdata.iteritems(),key= lambda x: int(x[1][0]), reverse=True):
        if data[4] == "NO_HIT":
            nhlengths.append(int(data[0]))
            nhgcs.append(float(data[1])/100)
            nhcoverages.append(float(data[2]))
        else:
            lengths.append(int(data[0]))
            gcs.append(float(data[1])/100)
            coverages.append(float(data[2]))
            orgindex = organisms.index(data[4])%num_colors
            color = sm.to_rgba(orgindex, alpha=float(data[3])/100)
            colors.append(color)                
            if data[4] in org_count:
                org_count[data[4]] += 1
            else:
                org_count[data[4]] = 1

##### METHOD TO CREATE RECTANGULAR PLOT
    xdata = np.array(gcs)*100
    ydata = np.array(coverages)
    area = np.array(lengths) #Double size
    nhxdata = np.array(nhgcs)*100
    nhydata = np.array(nhcoverages)
    nharea = np.array(nhlengths) #Double size


    #print "Total Contigs", len(lengths), len(nhlengths)

    #print r
    details_output.write("Minimum Coverage {}\n".format(ydata.min()))
    details_output.write("Mean Coverage {}\n".format(ydata.mean()))
    details_output.write("Maximum Coverage {}\n".format(ydata.max()))
    
    #print theta
    details_output.write("Minimum GC {}\n".format(xdata.min()))
    details_output.write("Mean GC {}\n".format(xdata.mean()))
    details_output.write("Maximum GC {}\n".format(xdata.max()))

    #print area
    details_output.write("Minimum Contig Size {}\n".format(area.min()))
    details_output.write("Mean Contig Size {}\n".format(area.mean()))
    details_output.write("Maximum Contig Size {}\n".format(area.max()))

    #Scaling Area and setting minimum size
    area = np.array(lengths)/50
    nharea = np.array(nhlengths)/50
    area = where( area < 10,10,area)
    nharea = where( nharea < 10,10,nharea)
    
    print "Creating figure"
    fig = figure(figsize=(11,11))

    details_output.write("Creating subplot\n")
    sp = fig.add_subplot(111,axisbg='none',frameon=True)
    sp.set_title("Contig GC,Length, Coverage and Taxonomy", fontsize=20, position=(0.5,1.0))
    sp.set_position([0.15,0.2,0.80,0.75])    

    sp.set_ylabel("Coverage")
    sp.set_xlabel("GC%")
    sp.set_xlim(0,100)
    sp.grid()

    if (not options.plotall):
        #Determine y-axis scaling (ignoring real outliers)
        pctl = 95
        while pctl < 99.9:
            #Why on earth can't I have this be pctl < 100.  For some reason it goes through the loop when int(pctl) == 100.
            #print "PCTL", pctl
            valueA = stats.scoreatpercentile(ydata, pctl)
            valueB = stats.scoreatpercentile(ydata, pctl+0.1)
            increase = (valueB-valueA)/valueB
            #print valueA, valueB, increase
            if (increase > 0.25):
                break
            pctl = pctl + 0.1
            #print "PCTL END", pctl
        max_y_value = stats.scoreatpercentile(ydata, pctl)
        details_output.write("Scaling maximum to {} percentile: {}\n".format(pctl,stats.scoreatpercentile(ydata, pctl)))
        sp.set_ylim(0,stats.scoreatpercentile(ydata, pctl))


    for index, cov in enumerate(ydata):
        if (cov > max_y_value ):
            ydata[index] = max_y_value
            details_output.write( "WARNING: Unit with coverage of {} printed at {} instead.  Out of scaling bounds.\n".format(cov,max_y_value))

    for index, cov in enumerate(nhydata):
        if (cov > max_y_value ):
            nhydata[index] = max_y_value
            details_output.write( "WARNING: Unit with coverage of {} printed at {} instead.  Out of scaling bounds.\n".format(cov,max_y_value))

    if (nhgcs):
        sp.scatter(nhxdata,nhydata,s=nharea,facecolor='none',edgecolor='gray',alpha=0.5,lw=1)    

    c = sp.scatter(xdata,ydata,c=colors,s=area,edgecolor='gray')
    #c = sp.scatter(xdata,ydata,c=colors,s=area,cmap=cmap, alpha=0.75,edgecolor='w')

    details_output.write("Setting colors\n")
    #sm = cm.ScalarMappable()
    #sm.set_cmap(cmap)
    #sm.set_array(c.get_array())
    #sm.autoscale()

    lpoints = []
    llabels = []

    p = Rectangle((0, 0), 1, 1, ec='gray', fc='none')
    lpoints.append(p)	
    llabels.append("No hit ["+str(len(nhlengths))+"]")

    for org in sorted(org_count,key= lambda org: org_count[org], reverse=True):
        if (len(llabels) == 30):
            details_output.write("WARNING: Over 30 different organisms.  Only showing top 30 in legend\n")
            break
        index = organisms.index(org)%num_colors
        color = sm.to_rgba(index)
        p = Rectangle((0, 0), 1, 1, fc=color)
        lpoints.append(p)
        llabels.append(org+" ["+str(org_count[org])+"]")
            

    fp = FontProperties(size='medium')
    sp.legend(lpoints,llabels,loc=9,frameon=1,bbox_to_anchor=(0.,-1.06,1.,1.),labelspacing=0.01,ncol=3,markerscale=0.5,prop=fp)
    figure_name = options.output + ".blast_bubbles.png"
    details_output.write("Creating {}\n".format(figure_name))
    details_output.close()
    savefig(figure_name)

    return 0

if __name__ == "__main__":
    sys.exit(main())
