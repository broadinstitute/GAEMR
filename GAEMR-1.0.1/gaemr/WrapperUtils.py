
# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.


import os
import shutil
import sys
import re
from multiprocessing import Process,Queue,Manager

import gaemr.PlatformConstant as pc
from gaemr.RunCommand import RunCommand
from gaemr.Site import Site
from gaemr.BamFile import BamFile

constant = pc.PlatformConstant()

#
# GAEMR wrapper utilities
#

class WrapperUtils(object):
    """This class represents components to make custom GAEMR wrappers."""

    def __init__(self, dir=None, force=False, output_header='assembly', threads = 4, extend=False, gaemr_dir_name='gaemr'):
        (self.gaemr_dir, self.work_dir, self.table_dir, self.chart_dir) = self.__setup_dirs(dir, force, extend, gaemr_dir_name)
        self.output_header = output_header
        self.pid_queue = Queue()
        self.threads = threads

    def get_gaemr_dir(self):
        return self.gaemr_dir

    def get_work_dir(self):
        return self.work_dir

    def get_chart_dir(self):
        return self.chart_dir

    def get_table_dir(self):
        return self.table_dir

    def get_output_header(self):
        return self.output_header

    def get_num_threads(self):
        return self.threads

    def get_pid_queue(self):
        return self.pid_queue

    def __setup_dirs(self, dir=None, force=False, extend=False, gaemr_dir_name=None):
        gaemr_dir = os.path.join(os.getcwd(),gaemr_dir_name)
        if dir:
            gaemr_dir = os.path.join(dir, gaemr_dir_name)

        if os.path.exists(gaemr_dir):
            if force:
                shutil.rmtree(gaemr_dir)
            elif not extend:
                print "GAEMR dir, " + str(gaemr_dir) + " exists.  Use --force to remove dir and run GAEMR analysis"
                sys.exit(-1)

        work_dir = os.path.join(gaemr_dir,'work')
        table_dir = os.path.join(gaemr_dir,'table')
        chart_dir = os.path.join(gaemr_dir,'chart')
        if not extend or not os.path.exists(gaemr_dir):
            os.mkdir(gaemr_dir)
            os.mkdir(work_dir)
            os.mkdir(table_dir)
            os.mkdir(chart_dir)

        return gaemr_dir,work_dir,table_dir,chart_dir

    def __concatenate_file_name(self, base, suffix):
        return base + '.' + suffix
    
    def __link_copy_files_helper(self, link_dir, arg_list, want_link=False):
        link_list = []

        for i in arg_list:
            link_name = None

            if i:
                link_name = os.path.join(link_dir,os.path.basename(i))
                if want_link:
                    os.symlink(os.path.abspath(i),link_name)
                else:
                    shutil.copyfile(os.path.abspath(i),link_name)
        
            link_list.append(link_name)
        
        return link_list

    def link_files(self, link_dir, arg_list):
        return self.__link_copy_files_helper(link_dir, arg_list, want_link=True)

    def copy_files(self, link_dir, arg_list):
        return self.__link_copy_files_helper(link_dir, arg_list, want_link=False)

    def print_start_of_process_pid(self):
        pid = str(os.getpid())
        print "[GAEMR] STARTED PID:  " + pid
        self.pid_queue.put(pid)

    def print_end_of_process_pid(self):
        pid = str(os.getpid())
        print "[GAEMR] COMPLETED PID:  " + pid
        self.pid_queue.put(pid)

    def print_and_run_command(self, rc):
        print "[GAEMR] PID:  " + str(os.getpid()) + " COMMAND:  " + rc.get_command()
        sys.stdout.flush()
        rc.run_command()

    def check_pids(self):
        pids = {}
        return_code = 0
        self.get_pid_queue().put('STOP')

        for pid in iter(self.get_pid_queue().get, 'STOP'):
            if pid not in pids:
                pids[pid] = 1
            else:
                pids[pid] += 1

        for pid in sorted(pids.iterkeys()):
            if int(pids[pid]) % 2:
                print "[GAEMR] FAILED PID:  " + pid
                return_code = -1
            else:
                print "[GAEMR] SUCCESSFUL PID:  " + pid

        return return_code

    def generate_html(self, assembly_name=None, html_subdir = 'html', base_url=None):
        if not assembly_name:
            assembly_name = self.get_output_header()
        Site(assembly_name, self.get_gaemr_dir(), html_subdir, base_url).build()

    def get_read_files(self, read_list_file=None):

        if read_list_file:
            read_file_dict = {}
            
            try:
                f = open(read_list_file,'r')
            except:
                print "Can't open read list file, " + read_list_file
                return -1

            for line in f:
                if re.match('#',line) or re.match('^$',line):
                    continue
                fields = line.rstrip('\n').split(',')
                if not fields[0] in read_file_dict:
                    read_file_dict[fields[0]] = {}
                if not fields[1] in read_file_dict[fields[0]]:
                    read_file_dict[fields[0]][fields[1]] = {}
                read_file_dict[fields[0]][fields[1]]['length'] = fields[2]
                read_file_dict[fields[0]][fields[1]]['direction'] = fields[3]
                read_file_dict[fields[0]][fields[1]]['insert_size'] = fields[4]
                if 'files' not in read_file_dict[fields[0]][fields[1]]:
                    read_file_dict[fields[0]][fields[1]]['files'] = []
                read_file_dict[fields[0]][fields[1]]['files'] += fields[5:]

            return read_file_dict            

    def check_read_files(self, read_file_dict=None):
        if read_file_dict:
            to_align_dict = {}

            for group in read_file_dict:
                to_align_dict[group] = {}
                for type in read_file_dict[group]:
                    to_align_dict[group][type] = {}
                    direction = read_file_dict[group][type]['direction']
                    to_align_dict[group][type]['file'] = self.revert_to_bam(read_file_dict[group][type]['files'],
                                                                            group, direction)
                    to_align_dict[group][type]['length'] = read_file_dict[group][type]['length']
                    to_align_dict[group][type]['insert_size'] = read_file_dict[group][type]['insert_size']
                    
            return to_align_dict

    def standardize_file_inputs(self, scaffolds=None, contigs=None, agp=None, minScaffSize=1, minConSize=1, minGapSize=10):
        if contigs and agp and scaffolds:
            return self.copy_files(self.get_work_dir(),[scaffolds,contigs,agp])

        self.print_start_of_process_pid()

        out_header = os.path.join(self.get_work_dir(),self.get_output_header())
        cmd_list = ['make_standard_assembly_files.py','-o',out_header,'--rename', '-s', str(minScaffSize),'-c', str(minConSize), '-g', str(minGapSize)]
        
        if contigs and agp:
            cmd_list += ['-C',contigs,'-A',agp]
        elif scaffolds:
            cmd_list += ['-S',scaffolds]
        elif contigs:
            cmd_list += ['-S',contigs]
        else:
            print "Must give a fasta file for make_standard_assembly_files."
            sys.exit(-1)

        rc = RunCommand(cmd_list)
        self.print_and_run_command(rc)
        self.print_end_of_process_pid()

        return self.__concatenate_file_name(out_header,"scaffolds.fasta"),\
               self.__concatenate_file_name(out_header, "contigs.fasta"), \
               self.__concatenate_file_name(out_header, "agp")

    def get_basic_assembly_stats(self, name=None, contigs=None, agp=None, assembler='assembler', extension='cumulative_sizes'):
        self.print_start_of_process_pid()

        if contigs and agp:
            if not name:
                name = self.get_output_header()
            output = os.path.join(self.get_table_dir(), name)
            chart = self.__concatenate_file_name(os.path.join(self.get_chart_dir(), name), extension)

            cmd_list = ['basic_assembly_stats.py','-n',name,'-a',assembler,'-o', output,'-C','-S','-t',chart,'-f', agp, contigs]
        
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)

            self.print_end_of_process_pid()

    def analyze_gap_ends(self, name=None, contigs=None, agp=None, extension='analyze_gap_ends'):
        self.print_start_of_process_pid()
        
        if contigs and agp:
            if not name:
                name = self.get_output_header()
            output = os.path.join(self.get_table_dir(), name)
            chart = self.__concatenate_file_name(os.path.join(self.get_chart_dir(), name), extension)

            cmd_list = ['analyze_gap_ends.py', '-c', chart, '-t', output, contigs, agp]

            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()

    def run_blast(self, query = None, db=constant.BLAST_NT, threads=None, extension='blast.xml', name=None, blast_task='megablast'):
        self.print_start_of_process_pid()

        if query:
            if not threads:
                threads=self.get_num_threads()
            if not name:
                name = self.get_output_header()
            blast_output = self.__concatenate_file_name(os.path.join(self.get_work_dir(), name),extension)
            
            rc = RunCommand(['run_blast.py','-o',blast_output,'-b',blast_task,'-t',str(threads), db, query])
            self.print_and_run_command(rc)

            self.print_end_of_process_pid()

            return blast_output

    def parse_blast_xml(self, blast_xml=None):
        self.print_start_of_process_pid()

        if blast_xml:
            parse_output = re.sub("xml", "parsed.txt", blast_xml)
            rc = RunCommand(['parse_blast_xml.py', '-o', parse_output, blast_xml])
            self.print_and_run_command(rc)
            
            self.print_end_of_process_pid()
            return parse_output

    def get_blast_hit_taxonomy(self, parsed_blast=None, query=None, nodes=constant.BLAST_NODES, names=constant.BLAST_NAMES, name=None):
        self.print_start_of_process_pid()

        if parsed_blast and query:
            if not name:
                name = self.get_output_header()
            tax_output = os.path.join(self.get_table_dir(), name)
            taxonomy_heatmap = re.sub("txt", "heatmap", parsed_blast)
            if nodes and names:
                rc = RunCommand(['get_blast_hit_taxonomy.py', '-o', tax_output,'-m', taxonomy_heatmap, parsed_blast, query])
                self.print_and_run_command(rc)

                self.print_end_of_process_pid()        

            return taxonomy_heatmap

    def blast_map(self, taxonomy_heatmap=None, agp=None, name=None):
        self.print_start_of_process_pid()

        if taxonomy_heatmap:
            if not name:
                name = self.get_output_header()
            blast_map_output = os.path.join(self.get_chart_dir(), name)

            cmd = ['blast_map.py','-o', blast_map_output]
            if agp:
                cmd += ['-g', agp]

	    cmd += [taxonomy_heatmap]
            rc = RunCommand(cmd)
            self.print_and_run_command(rc)
                
            self.print_end_of_process_pid()

    def analyze_assembly_rna(self, query=None, name=None):
        self.print_start_of_process_pid()

        if query:
            if not name:
                name = self.get_output_header()
            rnammer_output = os.path.join(self.get_work_dir(),self.__concatenate_file_name(name, 'rnammer.out'))
            rdp_output = os.path.join(self.get_work_dir(),self.__concatenate_file_name(name, 'rdp.out'))
            table_output = os.path.join(self.get_table_dir(),name)

            cmd_list = ['run_rna_analysis.py','-c','-f',rnammer_output,'-r',rdp_output,query]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)


            cmd_list = ['analyze_rna_hits.py','-c','-r',rdp_output,'-o',table_output,rnammer_output]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)

            self.print_end_of_process_pid()

    def run_nucmer(self, query=None, ref=None , name=None, extension='ref_vs_assembly'):
        self.print_start_of_process_pid()

        if query and ref:
            if not name:
                name = self.get_output_header()
            prefix = os.path.join(self.get_work_dir(), self.__concatenate_file_name(name, extension))
            cmd_list = ['run_nucmer.py','--mummerplot','-p',prefix, ref, query]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)

            self.link_files(self.get_chart_dir(), [self.__concatenate_file_name(prefix, 'png')])

            self.print_end_of_process_pid()
            
            return self.__concatenate_file_name(prefix,'coords')

    def compare_to_reference(self, coords_file=None, name=None):
        self.print_start_of_process_pid()

        if coords_file:
            if not name:
                name = self.get_output_header()
            table = os.path.join(self.get_table_dir(), name)
            cmd_list = ['compare_to_reference.py', '-o',table, '-n', '-c', coords_file]

            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            
            self.print_end_of_process_pid()

    def get_kmer_copy_number(self, fasta=None, name=None, kmer_size=29):
        self.print_start_of_process_pid()

        if fasta:
            if not name:
                name = self.get_output_header()
            cmd_list = ['kmer_copy_number.py', '-k', str(kmer_size), '-o', os.path.join(self.get_table_dir(), name), fasta]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
        
            self.print_end_of_process_pid()
        
    def run_kmer_coverage(self, ref=None, query=None, name=None, kmer_size=29):
        self.print_start_of_process_pid()

        if ref and query:
            if not name:
                name = self.get_output_header()
            cmd_list = ['kmer_coverage.py', '-k', str(kmer_size), '-o', os.path.join(self.get_table_dir(), name), ref, query]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            
            self.print_end_of_process_pid()
        
    def run_scaffold_accuracy(self, ref=None, query=None, name=None):
        self.print_start_of_process_pid()

        if ref and query:
            if not name:
                name = self.get_output_header()
            cmd_list = ['run_scaffold_accuracy.py', '-o', os.path.join(self.get_work_dir(), name), '-t', os.path.join(self.get_table_dir(), name),
                        ref, query]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
        
            self.print_end_of_process_pid()

    def __get_file_extension(self, file):
        return re.sub(".*\.","",file).upper()

    def __is_bam(self, file):
        return self.__get_file_extension(file) == 'BAM'

    def __is_aligned_bam(self, file):
        if self.__is_bam(file):
            return BamFile(file, self.__get_file_extension(file)).is_mapped()
        return 0

    def revert_to_bam(self, files=[], output_base='reads', direction='fr'):
        self.print_start_of_process_pid()

        if files:
            output_bam = os.path.join(self.get_work_dir(), self.__concatenate_file_name(output_base, 'unmapped.bam'))
            for i in files:
                if self.__is_bam(i) and not self.__is_aligned_bam(i):
                    os.symlink(i, output_bam)
                    self.print_end_of_process_pid()
                    return output_bam

            cmd_list = ['read_format_converter.py','-o',output_bam]
            if direction:
                cmd_list += ['-d', direction]
            cmd_list += files
                
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)

            self.print_end_of_process_pid()

            if not os.path.exists(output_bam):
                return None
            return output_bam

    def align_reads(self, unmapped_bam=None, ref=None, threads=None, aligner='BWA', ref_header='reference', make_index=True, align_type='-s'):
        self.print_start_of_process_pid()
        if unmapped_bam and ref:
            if not threads:
                threads = self.get_num_threads()
            output_header = re.sub("unmapped.bam", ref_header, unmapped_bam)
            cmd_list = ['align_reads.py', '-i', unmapped_bam, '-o', output_header, '-r', ref, '-a', aligner, align_type, '-t', str(threads), '-T', self.get_work_dir()]
            if not make_index:
                cmd_list += ['-x']
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)            
            self.print_end_of_process_pid()            

            return self.__concatenate_file_name(output_header, 'bam')

    def run_insert_size(self, bam_file=None, insert_size=None, std_dev=None):
        self.print_start_of_process_pid()

        if bam_file:
            output_header = re.sub(".bam", "", bam_file)
            cmd_list = ['run_insert_size_from_bam.py','-o',output_header]
            if insert_size:
                cmd_list += ['-i', str(insert_size)]
                if std_dev:
                    cmd_list += ['-s', str(std_dev)]
                   
                cmd_list += [bam_file]
                    
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)            
            self.print_end_of_process_pid()
            
            return self.__concatenate_file_name(output_header, 'insert_size.metrics')

    def plot_insert_size(self, insert_size_files=[], direction='fr', ref_header='reference', output_base='reads'):
        self.print_start_of_process_pid()

        if insert_size_files:
            type_output_header = os.path.join(self.get_table_dir(),self.__concatenate_file_name(output_base, ref_header))
            plot_output_header = os.path.join(self.get_chart_dir(),self.__concatenate_file_name(output_base, ref_header))

            cmd_list = ['plot_insert_size.py', '-o', plot_output_header,'-m', type_output_header, '-d', direction] + insert_size_files
            
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()

    def get_simple_bam_stats(self, bam_files=[], name=None, ref_header='reference'):
        self.print_start_of_process_pid()

        if bam_files:
            if not name:
                name = self.get_output_header()
            cmd_list = ['get_simple_bam_stats.py','-o', os.path.join(self.get_table_dir(),self.__concatenate_file_name(name, ref_header))] + bam_files
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()
            
    def get_bam_coverage_stats(self, bam_files=[], name=None, ref_header='reference', want_phys_cvg=False):
        self.print_start_of_process_pid()

        if bam_files:
            if not name:
                name = self.get_output_header()
            cmd_list = ['get_bam_coverage_stats.py','-o', os.path.join(self.get_table_dir(),self.__concatenate_file_name(name, ref_header))]
            if want_phys_cvg:
                cmd_list += ['-p']
            cmd_list += bam_files
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()


    def generate_bam_plots(self, bam_files=[], ref=None, name=None, ref_header='reference', window_size=1000):
        self.print_start_of_process_pid()

        if bam_files and ref:
            if not name:
                name = self.get_output_header()
            analysis_output_header = self.__concatenate_file_name(name, ref_header)
            data_dump_file = os.path.join(self.get_work_dir(),self.__concatenate_file_name(analysis_output_header, 'gc_cvg.details.txt'))
            plot_output_header = os.path.join(self.get_chart_dir(),analysis_output_header)
            histo_plot = self.__concatenate_file_name(plot_output_header,'gc_cvg')
            
            cmd_list = ["generate_bam_plots.py", "-g", ref, "-d", data_dump_file, "-w", str(window_size),
                        "-o", plot_output_header, "-hi", histo_plot] + bam_files

            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()

    def identify_coverage_anomalies(self, bam_files=[], name=None, ref_header='reference', window_size=1000):
        self.print_start_of_process_pid()

        if bam_files and ref_header:
            if not name:
                name = self.get_output_header()
            analysis_output_header = self.__concatenate_file_name(name, ref_header)
            coverage_anomalies = os.path.join(self.get_chart_dir(), self.__concatenate_file_name(analysis_output_header, "coverage_anomalies"))
            cmd_list = ['identify_coverage_anomalies.py','--window_size', str(window_size)] + bam_files + [coverage_anomalies]

            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()

    def make_detailed_table(self, name=None, aligned_bam_dict=None, contigs=None, agp=None, taxonomy_output=None):
        self.print_start_of_process_pid()

        # aligned_bam_dict=["type" =
        #                    ["group1 = [
        #                       [file1, file2, filen]
        #                                           
        #                                ]
        #                    ]
        #                 ]
        if contigs and agp:
            if not name:
                name = self.get_output_header()
            out_header = os.path.join(self.get_table_dir(),name)
            cmd_list = ['make_detailed_assembly_table.py','-c',out_header,'-s',out_header,'-a',agp]
            if taxonomy_output:
                cmd_list += ['-t',taxonomy_output]

            if aligned_bam_dict:
                for type in aligned_bam_dict.keys():
                    arg = "--" + str(type) + "_bam"
                    bam_file_string = ''

                    for group in aligned_bam_dict[type].keys():
                        bam_file_string += aligned_bam_dict[type][group]['file'] + ','
                    cmd_list += [arg,bam_file_string[:-1]]

            cmd_list += [contigs]

            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            self.print_end_of_process_pid()
            return self.__concatenate_file_name(out_header, 'contig_detail.table.txt')
        
    def blast_bubbles(self, name=None, taxonomy_output=None, contig_detail=None):
        self.print_start_of_process_pid()
        
    
        if taxonomy_output and contig_detail:
            if not name:
                name = self.get_output_header()
            blast_bubble_output = os.path.join(self.get_chart_dir(),name)
            blast_bubble_detail = os.path.join(self.get_work_dir(),name)
            cmd_list = ['blast_bubbles.py', '-v', blast_bubble_detail,'-o', blast_bubble_output, contig_detail, taxonomy_output]
            rc = RunCommand(cmd_list)
            self.print_and_run_command(rc)
            
            self.print_end_of_process_pid()
