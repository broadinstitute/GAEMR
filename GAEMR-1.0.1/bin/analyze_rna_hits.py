#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import sys
from optparse import OptionParser
from gaemr.RnaClass import RnammerFile,RDPFile
import gaemr.PlatformConstant as pc
constant = pc.PlatformConstant()
from gaemr.SimpleTable import SimpleTable
    
parser = OptionParser(usage="usage: %prog [options] rnammer output fasta")

parser.add_option('--classify', '-c', action="store_true", default=None,dest="classify",
                  help='Whether or not to print taxonomy of rna hits (default=%default)')
parser.add_option('--rdp_classifier_output', '-r', action="store", default=None, type='string',dest="rdp_out",
                  help='RDP Classifier output filename (default=%default)')
parser.add_option('--gene', '-g', action="store", default="all", type='string', dest="gene",
                  help='Gene of interest: 5s, 16s, 23s or all (default=%default)')
parser.add_option('--lineage', '-l', action="store", default="genus", type='string', dest="lineage",
                  help='Level of lineage to print out in classification output (default=%default)')
parser.add_option('--score', '-s', action="store", default="0.0", type='float', dest="score",
                  help='Cutoff score for reporting taxonomic hits in classification output - values from 0 to 1 (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
                  help='Output header for results (default=%default)')
parser.add_option('--no_html',action="store_false",default=True,dest="html",
                  help='Whether or not to print out html table output (default=%default)')

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.error("Must supply rnammer fasta file as input.")
    sys.exit(1)

# private function to make sure user put in right args
def __check_inputs(classify,rdp_out):
    if classify and not rdp_out or rdp_out and not classify:
        return 0
    return 1

def __summarize_results(s_dict,g_count):
    title = 'rRNA Analysis Summary'
    headers = ['Gene','Total Copies','Lineage','Number Organisms Found','Organism IDs']
    data = []
    for gene in s_dict:
        orgs = sorted(s_dict[gene].keys())
        string = ''.join([orgs[i] + ';' for i in xrange(len(orgs))])
        data.append([gene, g_count[gene],options.lineage, len(orgs), string[:-1]])

    return headers,data,title

def __print_rna_analysis(molecules):
    # print out the info
    title = 'rRNA Analysis'
    rna_table = []
    rna_table_headers = ['Gene','Genus','Query ID','Query Start','Query End','Hit Direction','Length']
    if options.classify:
        rna_table_headers.append('Taxonomy')

    rna_count = 0

    summary_dict = {}
    gene_count = {}
    for m in sorted(molecules.iterkeys()):
        for i in molecules[m]:
            rna_table.append([])
            gene = i.get_molecule()
            lineage = i.get_genus()
            if options.gene == "all":
                rna_table[rna_count] += [gene, i.get_genus(), i.get_contig(), i.get_start(),
                                         i.get_end(), i.get_direction(), i.get_length()]
                if options.classify:
                    if gene not in summary_dict:
                        summary_dict[gene] = {}
                        gene_count[gene] = 0
                    if lineage not in summary_dict[gene]:
                        summary_dict[gene][lineage] = 0
                    summary_dict[gene][lineage] += 1
                    gene_count[gene] += 1

            else:
                if m == options.gene:
                    rna_table[rna_count] += [i.get_molecule(), i.get_genus(), i.get_contig(), i.get_start(),
                                             i.get_end(), i.get_direction(), i.get_length()]
                if options.classify:
                    if gene not in summary_dict:
                        summary_dict[gene] = {}
                        gene_count[gene] = 0
                    if lineage not in summary_dict[gene]:
                        summary_dict[gene][lineage] = 0
                    summary_dict[gene][lineage] += 1
                    gene_count[gene] += 1

            if options.classify:
                rna_table[rna_count] += [i.get_taxonomy()]

            rna_count += 1

    st = SimpleTable(rna_table_headers,rna_table,title)
    output = None
    if options.output:
        output = options.output + ".rna_analysis_details"
    st.print_output(output,options.html)



    if options.classify:
        headers, data, title = __summarize_results(summary_dict,gene_count)
        output = None
        if options.output:
            output = options.output + ".rna_analysis_summary"
        st = SimpleTable(headers,data,title)
        st.print_output(output,options.html)

def main():
    if __check_inputs(options.classify, options.rdp_out):

        rf_obj = RnammerFile(args[0])
        molecules = rf_obj.get_molecules()

        if options.classify:
            rdp_obj = RDPFile(options.rdp_out, molecules, options.lineage, options.score)

        __print_rna_analysis(molecules)

        return 0
    else:
        print "If classifying hits, must supply classify flag and rdp output file."
        sys.exit(-1)

if __name__ == "__main__":
    sys.exit(main())
