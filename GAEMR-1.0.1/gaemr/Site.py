#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import os.path
import os
import glob
import shutil
import time


#########################################################
# this class generates HTML for GAEMR                   #
# output. It is tightly bound to the                    #
# current (as of 8/17) output files of                  #
# gaemr and any modification to that                    #
# required modifications here.                          #
#                                                       #
#This requires 3 inputs                                 #
#                                                       #
# 1. An assembly name (for prefixes)                    #
# 2. The location of GAEMR output                       #
# 3. Where you want to build the site (subdir of gaemr) #
#                                                       #
#To use it, 2 steps are needed.                         #
# 1. First build object                                 #
#   for example:                                        #
#   site = Site('G15640_allpaths_200x100x0x_3009',      #
#               '/Users/harindra/Desktop/gaemr',        #
#               'html')                                 #
#                                                       #
# 2. Call "build" method to build                       #
#   site.build()                                        #
#########################################################


#generate HTML site for GAEMR output
class Site:
    """generate HTML site for GAEMR output"""

    #start here
    def __init__(self,assembly_name,gaemr_dir,site_subdir='html',home_base=None):
        """initialize class"""
        self.gaemr_dir = gaemr_dir
        self.assembly_name = assembly_name
        self.site_dir=os.path.join(gaemr_dir, site_subdir)
        self.template_dir=os.path.join(os.path.dirname(__file__),'resources/templates').replace('\\','/')
        self.site_image_dir = os.path.join(os.path.dirname(__file__),'resources/img').replace('\\','/')
        self.site_css_dir = os.path.join(os.path.dirname(__file__),'resources/css').replace('\\','/')
        self.create_file_system_structure_for_site()
        self.home_base = home_base

        os.chdir(self.site_dir)

    #software version
    def get_version(self):
        """get software version"""
        return "1.8_10.24.2012"

    #build the site now
    def build(self):
        """build the site now,"""
        print "generating assembly stats section ...."
        html = self.get_section_assembly_stats()
        self.write_html_to_file(html,'assembly_stats.html')

        print "generating read alignment stats section ...."
        html = self.get_section_read_alignment_stats()
        self.write_html_to_file(html,'read_stats.html')

        print "generating Kmer analysis sections ...."
        html = self.get_section_kmer_analysis()
        self.write_html_to_file(html,'kmer_analysis.html')

        print "generating Gap end analysis sections ...."
        html = self.get_section_gap_end_analysis()
        self.write_html_to_file(html,'gap_end_analysis_summary.html')

        html = self.get_section_reference_analysis()
        if html:
            print "generating reference analysis section ...."
            self.write_html_to_file(html,'reference_analysis.html')
        else:
            print "skipping reference analysis section ...."

        print "generating blast analysis section ...."
        html = self.get_section_blast_analysis()
        self.write_html_to_file(html,'blast_analysis.html')

        print "generating rRNA Details section ...."
        html = self.get_rrna_details()
        self.write_html_to_file(html,'rrna_details.html')

        print "generating Assembly Details section ...."
        html = self.get_assembly_details()
        self.write_html_to_file(html,'assembly_details.html')

        print "generating gaemr summary page ...."
        html = self.get_gaemr_summary()
        self.write_html_to_file(html, 'gaemr_summary.html')
        
        print "generating main index page ...."
        html = self.get_main_index()
        self.write_html_to_file(html,'index.html')



    #fetches HTML output from gaemr module
    def __fetch_html_from_file__(self,file_name):
        """fetches HTML output from gaemr module"""
        err_html = "<h3>sorry,that information was not generated in this run </h3>"
        content=''
        if file_name=='-1':
            return err_html
        else:
            html_content = open(file_name,"r")
            for line in html_content:
                content += line
            html_content.close()
        if len(content)==0:
            return "<h3>sorry, "+file_name + ' was empty in this run </h3>'
        return content


    def get_gaemr_summary(self):
        """generate HTML for gaemr summary section"""
        
        #read template in
        template_file = open(self.template_dir+'/gaemr_summary.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()        

        print "\t - generating figure cumulative sizes image ...."
        cumul_sizes_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.cumulative_sizes.png')#        

        print "\t - generating blast bubbles image ...."
        blast_bubbles_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.blast_bubbles.png')

        print "\t - generating coverage line plot and GC% heat map combined plot file ...."
        scaffold_cvrg_line_and_gc_heat_images=self.__find_single_file_that_look_like_this__('../chart/' + '*.scaffolds.gc_cvg.png')

        print "\t - generating rRNA analysis summary ...."
        rna_analysis_sum_html_content_file_name = self.__find_single_file_that_look_like_this__('../table/'+self.assembly_name + '.rna_analysis_summary.table.html')
        rna_analysis_sum_html_content_to_inject =self.__fetch_html_from_file__(rna_analysis_sum_html_content_file_name)

        print "\t - generating html for ref vs assembly image ...."
        ref_vs_assembly_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.ref_vs_assembly.png')#        

        substitution_table = {
            'cumulative_sizes_image':cumul_sizes_image,
            'blast_bubbles_image':blast_bubbles_image,
            'scaffold_gc_heat_cvg_line_images':scaffold_cvrg_line_and_gc_heat_images,
            'rna_analysis_summary':rna_analysis_sum_html_content_to_inject,
            'ref_vs_assembly':ref_vs_assembly_image,
            }

        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)        


    #generate HTML for Section 1: Assembly Stats
    def get_section_assembly_stats(self):
        """generate HTML for assembly stats section"""

        #read template in
        template_file = open(self.template_dir+'/assembly_stats.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        print "\t - generating basic assembly stats table ...."
        basic_assembly_stat_html_content_to_inject =\
        self.__fetch_html_from_file__(
            '../table/'+self.assembly_name + '.basic_assembly_stats.table.html')

        print "\t - generating rRNA analysis summary ...."
        file_name_of_html = self.__find_single_file_that_look_like_this__('../table/'+self.assembly_name + '.rna_analysis_summary.table.html')
        rna_analysis_sum_html_content_to_inject = self.__fetch_html_from_file__(file_name_of_html)

        print "\t - generating figure cumulative sizes image ...."
        cumul_sizes_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.cumulative_sizes.png')

        substitution_table  = {
            'basic_assembly_stats':basic_assembly_stat_html_content_to_inject,
            'rna_analysis_summary':rna_analysis_sum_html_content_to_inject,
            'cumulative_sizes_image':cumul_sizes_image,
            }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)


    #find files that look like this in the FS
    def __find_files_that_look_like_this__(self,pattern):
        """look for files in the FS"""
        return glob.glob(pattern)


    #find single file that look like this in the FS
    def __find_single_file_that_look_like_this__(self,pattern):
        """look for files in the FS"""
        files = glob.glob(pattern)
        if len(files)>0:
            return files[0]
        else:
            print "\t\t----WARNING: no file (single file search) found for pattern: "+ pattern + " ----"
            return '-1'

    #if HTML content is missing in the data input, return a HTML
    #message saying so in the string, if content found, return as is
    def __add_missing_data_message__(self,found_data,content_name):
        """if HTML content is missing in the data input, return a HTML message saying so"""
        if len(found_data) == 0:
            print "\t\t----WARNING: no HTML found for "+ content_name + "----"
            return "<h3>sorry, "+content_name+ " information was not generated in this run </h3>"
        else:
            return found_data


    #django free - in progr
    #generate HTML for Section 2: Read alignment Stats
    def get_section_read_alignment_stats(self):
        """generate HTML for \"Section 2: Read alignment Stats\""""

        #read template in
        template_file = open(self.template_dir+'/read_stats.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        #dynamic content

        print "\t - generating simple BAM stats table ...."
        simple_bam_stat_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.scaffolds.simple_bam_stats.table.html'):
            simple_bam_stat_html_to_inject+=self.__fetch_html_from_file__(file)
        print "\t - generating simple BAM sequence coverage table ...."
        seq_cvrg_stats_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.scaffolds.bam_seq_cvg_stats.table.html'):
            seq_cvrg_stats_html_to_inject +=self.__fetch_html_from_file__(file)
        print "\t - generating simple BAM physical coverage table ...."
        phys_cvrg_stats_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.scaffolds.bam_phys_cvg_stats.table.html'):
            phys_cvrg_stats_html_to_inject +=self.__fetch_html_from_file__(file)
        print "\t - generating simple BAM insert size table ...."
        insert_size_stats_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.scaffolds.insert_size.table.html'):
            insert_size_stats_html_to_inject +=self.__fetch_html_from_file__(file)
        print "\t - generating scaffold insert size figure ...."
        scaf_insert_size_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.insert_size.png')
        print "\t - generating scaffold coverage anomalie figure ...."
        scaf_cov_anorm_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.coverage_anomalies_cvrg_anomalie_scatter.png')
        print "\t - generating scaffold coverage anomalie cluster length figure ...."
        scaf_cov_anorm_clust_len_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.coverage_anomalies_cvrg_anomalie_cluster_length_bar.png')
        print "\t - generating gc coverage figure ...."
        gc_cvg_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.gc_cvg*.histogram.png')
        print "\t - generating coverage change delta figure (scatter plot) ...."
        cv_change_scatter_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.coverage_anomalies.coverage_change_across_genome_scatter.png')
        print "\t - generating coverage change delta figure (histogram plot) ...."
        cv_change_histogram_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.coverage_anomalies.coverage_change_across_genome_histogram.png')

        print "\t - generating coverage line plot and GC% heat map combined plot file ...."
        scaffold_cvrg_line_and_gc_heat_images=self.__find_files_that_look_like_this__('../chart/' + '*.scaffolds.gc_cvg.png')


        substitution_table_html = {
            'simple_bam_stats':self.__add_missing_data_message__(simple_bam_stat_html_to_inject,'simple bam stat table'),
            'seq_cvrg_stats':self.__add_missing_data_message__(seq_cvrg_stats_html_to_inject,'sequence coverage table'),
            'phys_cvrg_stats':self.__add_missing_data_message__(phys_cvrg_stats_html_to_inject,'physical coverage table'),
            'insert_size_stats':self.__add_missing_data_message__(insert_size_stats_html_to_inject,'insert size table'),
                                   }
        substitution_table_images= {
            'scaf_inserts_images':scaf_insert_size_images,
            'scaf_cov_anorm':scaf_cov_anorm_images,
            'scaf_cov_anorm_clust_len_images':scaf_cov_anorm_clust_len_images,
            'gc_cvg_images':gc_cvg_images,
            'coverage_change_scatter_images':cv_change_scatter_images,
            'coverage_change_histogram_images':cv_change_histogram_images,
            'scaffold_gc_heat_cvg_line_images':scaffold_cvrg_line_and_gc_heat_images,
            }
        return self.add_this_muliple_value_substitution_info_table_to_html(template_content,substitution_table_html,substitution_table_images)



    #DJANGO FREE
    #get the HTML output of gaemr modules related to kmer analysis
    def get_section_kmer_analysis(self):
        """get the HTML output of gaemr modules related to kmer analysis"""

        #read template in
        template_file = open(self.template_dir+'/kmer_analysis.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        #dynamic content
        print "\t - generating kmer copy number table ...."
        file=self.__find_single_file_that_look_like_this__('../table/' + '*.kmer_copy_number.table.html')#
        kmer_copy_number_html_to_inject=self.__fetch_html_from_file__(file)
        print "\t - generating kmer coverage table ...."

        substitution_table = {
            'kmer_copy_number':self.__add_missing_data_message__(kmer_copy_number_html_to_inject,'kmer copy number table')
                            }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)

    #Django free
    #get the Gap End analysis Section 4
    def get_section_gap_end_analysis(self):
        """generate hTML for gap end analysis"""

        #read template in
        template_file = open(self.template_dir+'/gap_end_analysis.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        print "\t - generating gap end analysis table ...."
        table_file=self.__find_single_file_that_look_like_this__('../table/' + '*.gap_analysis.table.html')#
        gap_analysis_table_html_to_inject=self.__fetch_html_from_file__(table_file)

        print "\t - generating gap end simple sequence table ..."
        table_file=self.__find_single_file_that_look_like_this__('../table/' + '*.gap_ss_analysis.table.html')
        ss_analysis_table_html_to_inject=self.__fetch_html_from_file__(table_file)
        
        print "\t - generating gap end analysis PNG files ...."
        cg_sizes_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.cg_sizes.png')#
        cg_distinctness_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.cg_distinctness.png')#
        cg_copy_num_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.cg_copy_number.png')#
        cg_gc_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.cg_gc.png')#
        ue_distinctness_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.ue_distinctness.png')#
        ue_copy_num_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.ue_copy_number.png')#
        ue_gc_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.analyze_gap_ends.ue_gc.png')#
        substitution_table = {
            'gap_analysis_table':self.__add_missing_data_message__(gap_analysis_table_html_to_inject,'gap analysis table'),
            'gap_ss_table':self.__add_missing_data_message__(ss_analysis_table_html_to_inject,'gap simple sequence table'),
            'cg_sizes_image':cg_sizes_image,
            'cg_distinctness':cg_distinctness_image,
            'cg_copy_num':cg_copy_num_image,
            'cg_gc':cg_gc_image,
            'ue_distinctness':ue_distinctness_image,
            'ue_copy_num':ue_copy_num_image,
            'ue_gc':ue_gc_image,
            }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)

    #DJANGO FREE
    #generate HTML for blast analysis section
    def get_section_blast_analysis(self):
        """generate HTML for blast analysis section"""

        #read template in
        template_file = open(self.template_dir+'/blast_analysis.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        #find out dynamic  content names and files
        print "\t - generating blast taxonomy details table ...."
        table_file=self.__find_single_file_that_look_like_this__('../table/' + '*.blast_hit_taxonomy.table.html')
        blast_tax_det_table_html_to_inject=self.__fetch_html_from_file__(table_file)
        print "\t - generating blast bubbles image ...."
        blast_bubbles_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.blast_bubbles.png')
        print "\t - generating blast heat map image ...."
        blast_heatmap_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.blast_map.png')

        substitution_table = {
            'blast_tax_det_table':blast_tax_det_table_html_to_inject,
            'blast_bubbles_image':blast_bubbles_image,
            'blast_heatmap_image':blast_heatmap_image,
            }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)

    #DJANGO FREE
    #generate HTML for rRNA Details section
    def get_rrna_details(self):
        """generate HTML for rRNA Details section"""

        #read template in
        template_file = open(self.template_dir+'/rrna_details.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        #get specific file details
        print "\t - generating rRNA analysis details table ...."
        table_file=self.__find_single_file_that_look_like_this__('../table/' + '*.rna_analysis_details.table.html')
        rrna_analysis_details_html_to_inject=self.__fetch_html_from_file__(table_file)


        substitution_table = {
                                'rrna_analysis_details':rrna_analysis_details_html_to_inject
                             }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)

    #DJANGO FREE
    #generate HTML for assembly details section
    def get_assembly_details(self):
        """generate HTML for assembly details section"""

        #read template in
        template_file = open(self.template_dir+'/assembly_details.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()


        print "\t - generating contig details table ...."
        table_file=self.__find_single_file_that_look_like_this__('../table/' + '*.contig_detail.table.html')#
        contig_details_html_to_inject=self.__fetch_html_from_file__(table_file)

        print "\t - generating scaffold details table ...."
        table_file=self.__find_single_file_that_look_like_this__('../table/' + '*.scaffold_detail.table.html')#
        scaffold_details_html_to_inject=self.__fetch_html_from_file__(table_file)
        substitution_table = {
            'contig_details':contig_details_html_to_inject,
            'scaffold_details':scaffold_details_html_to_inject,
            }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)


    #get reference analysis section HTML
    def get_section_reference_analysis(self):
        """reference analysis section"""
        
        # use this as a canonical test for existance of reference analysis, and skip section
        if not self.__find_files_that_look_like_this__('../chart/' + '*.ref_vs_assembly.png'):
            return None

        #read template in
        template_file = open(self.template_dir+'/reference_analysis.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        #add dynamic content to template next, find files first

        print "\t - generating scaffold accuracy table ...."
        scaf_accuracy_table=self.__find_single_file_that_look_like_this__('../table/' + '*.scaffold_accuracy.table.html')#
        scaf_accuracy_table_html_to_inject=self.__fetch_html_from_file__(scaf_accuracy_table)

        print "\t - generating html for ref vs assembly image ...."
        ref_vs_assembly_image=self.__find_single_file_that_look_like_this__('../chart/' + '*.ref_vs_assembly.png')#

        file=""
        print "\t - generating reference simple bam stat files ...."
        ref_simple_bam_stat_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.reference.simple_bam_stats.table.html'):
            ref_simple_bam_stat_html_to_inject += self.__fetch_html_from_file__(file)
        ref_insert_size_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.insert_size.png')

        file=""
        print "\t - generating reference comparison To reference file...."
        compare_to_ref_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.compare_to_ref.table.html'):
            compare_to_ref_html_to_inject += self.__fetch_html_from_file__(file)

        print "\t - generating sequence coverage files ...."
        ref_seq_cvrg_stat_table_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.reference.bam_seq_cvg_stats.table.html'):
            ref_seq_cvrg_stat_table_html_to_inject += self.__fetch_html_from_file__(file)
        cvrg_anormalie_scatter_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.coverage_anomalies_cvrg_anomalie_scatter.png')

        print "\t - generating physical coverage files ...."
        phys_cvrg_stat_table_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.reference.bam_phys_cvg_stats.table.html'):
            phys_cvrg_stat_table_html_to_inject += self.__fetch_html_from_file__(file)
        cluster_length_bar_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.coverage_anomalies_cvrg_anomalie_cluster_length_bar.png')

        print "\t - generating coverage anomalie (change in coverage) scatter plot files ...."
        ref_change_in_cvrg__scatter_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.coverage_anomalies.coverage_change_across_genome_scatter.png')

        print "\t - generating coverage anomalie (change in coverage) histogram plot files ...."
        ref_change_in_cvrg__histogram_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.coverage_anomalies.coverage_change_across_genome_histogram.png')

        print "\t - generating insert size files ...."
        insert_size_table_html_to_inject=""
        for file in self.__find_files_that_look_like_this__('../table/' + '*.reference.insert_size.table.html'):
            insert_size_table_html_to_inject += self.__fetch_html_from_file__(file)
        insert_size_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.gc_cvg*.histogram.png')

        print "\t - generating kmer coverage file ...."
        file=self.__find_single_file_that_look_like_this__('../table/' + '*.kmer_coverage.table.html')#
        kmer_coverage_html_to_inject=self.__fetch_html_from_file__(file)

        print "\t - generating coverage line plot and GC% heat map combined plot file ...."
        ref_cvrg_line_and_gc_heat_images=self.__find_files_that_look_like_this__('../chart/' + '*.reference.gc_cvg.png')

        substitution_table_html = {
            'scaffold_accuracy_table':self.__add_missing_data_message__(scaf_accuracy_table_html_to_inject,'scaffold accuracy table'),
            'ref_simple_bam_stat':self.__add_missing_data_message__(ref_simple_bam_stat_html_to_inject,'reference simple bam stat'),
            'ref_seq_cvrg_stat':self.__add_missing_data_message__(ref_seq_cvrg_stat_table_html_to_inject,'ref sequence cvrg stat'),
            'phys_cvrg_stat':self.__add_missing_data_message__(phys_cvrg_stat_table_html_to_inject,'phys cvrg stat'),
            'insert_size_table':self.__add_missing_data_message__(insert_size_table_html_to_inject,'insert size table'),
            'kmer_coverage':self.__add_missing_data_message__(kmer_coverage_html_to_inject,'kmer coverage table'),
            'comparison_to_reference':self.__add_missing_data_message__(compare_to_ref_html_to_inject,'comparison to reference'),
            }

        #####
        # NOTE: substitution function that replaces template tag with below content expects the value to be a list,
        # hence even if value is only a single image,STILL enclose in list
        substitution_table_images = {
            'ref_vs_assembly':[ref_vs_assembly_image],
            'ref_insert_size':ref_insert_size_images,
            'cvrg_anormalie_scatter':cvrg_anormalie_scatter_images,
            'cluster_length_bar':cluster_length_bar_images,
            'ref_change_in_cvrg__scatter_images':ref_change_in_cvrg__scatter_images,
            'ref_change_in_cvrg__histogram_images':ref_change_in_cvrg__histogram_images,
            'insert_size_image':insert_size_images,
            'cvrg_and_gc_heatmap':ref_cvrg_line_and_gc_heat_images,
            }
        return self.add_this_muliple_value_substitution_info_table_to_html(template_content,substitution_table_html,substitution_table_images)

    #copy over helper files used by index.html
    def __copy_over_helper_html_files__(self):
        """copy over helper files used by index.html"""
        print "\t - copying over helper HTML file to assembly_analysis_ind.html" + self.site_dir + " ...."
        shutil.copyfile(self.template_dir+'/assembly_analysis_ind.html',self.site_dir+'/assembly_analysis_ind.html')

        print "\t - copying over helper HTML file to assembly_analysis_title.html" + self.site_dir + " ...."
        #shutil.copyfile(self.template_dir+'/assembly_analysis_title.html',self.site_dir+'/assembly_analysis_title.html')
        html = self.get_assembly_analysis_title_html()
        self.write_html_to_file(html,'assembly_analysis_title.html')

        print "\t - copying over broad logo to " + self.site_dir + "/img ...."
        os.makedirs(self.site_dir+'/img')
        shutil.copyfile(self.site_image_dir+'/broad-logo.jpg',self.site_dir+'/img/broad-logo.jpg')
        print "\t - copying over css to " + self.site_dir + "/css ...."
        os.makedirs(self.site_dir+'/css')
        shutil.copyfile(self.site_css_dir+'/gaemr.css',self.site_dir+'/css/gaemr.css')


    #Django free
    #generate the main index page
    def get_assembly_analysis_title_html(self):
        """generate the get_assembly_analysis_title page"""

        #read template in
        template_file = open(self.template_dir+'/assembly_analysis_title.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()

        title_assembly_name=''
        if self.home_base != None:
            title_assembly_name += '<a href="javascript:goTo(' + "'"+ self.home_base + "'" +');"/>' + self.assembly_name +'</a>'
        else:
            title_assembly_name += self.assembly_name
        #add dynamic content
        localtime = time.asctime( time.localtime(time.time()) )
        substitution_table = {'assembly_name': title_assembly_name,
                              'site_generation_time':localtime,
                                }
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)


    #Django free
    #generate the main index page
    def get_main_index(self):
        """generate the main index page"""
        self.__copy_over_helper_html_files__()

        #read template in
        template_file = open(self.template_dir+'/index.html','r')
        template_content=''
        for line in template_file:
            template_content += line
        template_file.close()
        substitution_table = {}
        return self.add_this_single_value_substitution_info_table_to_html(template_content,substitution_table)



    #create a HTML file using given HTML and file name in given dir
    def write_html_to_file(self,html,file_name):
        """create a HTML file using given HTML and file name in given dir"""
        html_write = open(self.site_dir+'/'+file_name,'w')
        html_write.write(html)
        html_write.close()

    #creates the file system structure for gaemr site
    def create_file_system_structure_for_site(self):
        """creates the file system structure for gaemr site"""
        if os.path.exists(self.site_dir):
            print "warning: an older gaemr site exists....deleting...."
            shutil.rmtree(self.site_dir)
        os.makedirs(self.site_dir)


    #generic function to substitute a keyword with key/value where value is always a single item
    def add_this_single_value_substitution_info_table_to_html(self,html_template,substitution_table):
        """ adds a single value to a HTML file """
        html=html_template
        for key,value in substitution_table.iteritems():
            #print "KEY:  ", key, " VALUE: ", value
            tag = '{{ '+key+' }}'
            html = html.replace(tag,value)
        return html




    #generic function to substitute a keyword with key/value, where value can be a list of tables (html) or images (png)
    def add_this_muliple_value_substitution_info_table_to_html(self,html_template,substitution_html_tables,substitution_images):
        """ adds a multi value list to a HTML file """
        html=html_template

        #replace all html tags in template with html content given in table
        for key,tables in substitution_html_tables.iteritems():
            tag = '{{ '+key+' }}'
            html_of_all_tables_in_this_key=''
            for table in tables:
                html_of_all_tables_in_this_key += table

            html = html.replace(tag,html_of_all_tables_in_this_key)

        #add the images
        for key,images in substitution_images.iteritems():
            tag = '{{ '+key+' }}'
            #generate all img html tags
            html_img_tags=''
            for image in images:
                img_tag = '<img src="' + image + '" alt="reference related image missing"'+'class="plot"/>'
                html_img_tags += img_tag + '\n'
            #add img tags to html template
            html = html.replace(tag,html_img_tags)

        return html
