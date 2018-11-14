#
# Helper class to handle the basic assembly stats
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import PlatformConstant as pc
constant = pc.PlatformConstant()

# Fixing the order of output of stats
C_STATS = ['Contigs','Max Contig','Mean Contig','Contig N50','Contig N90','Total Contig Length','Assembly GC']
S_STATS = ['Scaffolds','Max Scaffold','Mean Scaffold','Scaffold N50','Scaffold N90','Total Scaffold Length']
G_STATS = ['Captured Gaps','Max Gap','Mean Gap','Gap N50','Total Gap Length']

# data structures to hold the stats being reported above
CON_STATS = {'Contigs':[],
             'Max Contig':[],
             'Mean Contig':[],
             'Contig N50':[],
             'Contig N90':[],
             'Total Contig Length':[],
             'Assembly GC':[]
             }
SCAFF_STATS = {'Scaffolds':[],
               'Max Scaffold':[],
               'Mean Scaffold':[],
               'Scaffold N50':[],
               'Scaffold N90':[],
               'Total Scaffold Length':[]
               }
GAP_STATS = {'Captured Gaps':[],
             'Max Gap':[],
             'Mean Gap':[],
             'Gap N50':[],
             'Total Gap Length':[]
             }

# default delimiter
DELIMITER = constant.TABLE_DELIMITER

class AssemblyStatsUtil:
   """This class handles some of the assembly stats reporting."""

   def __init__(self,assembly):
       self.assembly = assembly
       num_assemblies = len(assembly)
       self.__load_contig_stats(assembly)
       self.__load_scaffold_stats(assembly)
       self.__load_gap_stats(assembly)


   # private contig stats loader
   def __load_contig_stats(self,assembly):
      """Populates the contig stats"""
      for a in assembly:
         CON_STATS['Contigs'].append(a.num_contigs())
         CON_STATS['Max Contig'].append(a.get_max_contig())
         CON_STATS['Mean Contig'].append(a.contig_mean())
         CON_STATS['Contig N50'].append(a.contigN50())
         CON_STATS['Contig N90'].append(a.contigN90())
         CON_STATS['Total Contig Length'].append(a.get_ungapped_length())
         CON_STATS['Assembly GC'].append(a.get_assembly_GC())

   # private scaffold stats loader
   def __load_scaffold_stats(self,assembly):
      """Populates the scaffold stats"""
      for a in assembly:
         SCAFF_STATS['Scaffolds'].append(a.num_scaffolds())
         SCAFF_STATS['Max Scaffold'].append(a.get_max_scaffold())
         SCAFF_STATS['Mean Scaffold'].append(a.scaffold_mean())
         SCAFF_STATS['Scaffold N50'].append(a.scaffoldN50())
         SCAFF_STATS['Scaffold N90'].append(a.scaffoldN90())
         SCAFF_STATS['Total Scaffold Length'].append(a.get_gapped_length())

   # private gap stats loader
   def __load_gap_stats(self,assembly):
      """Populates the gap stats"""
      for a in assembly:
         GAP_STATS['Captured Gaps'].append(a.num_gaps())
         GAP_STATS['Max Gap'].append(a.gap_max())
         GAP_STATS['Mean Gap'].append(a.gap_mean())
         GAP_STATS['Gap N50'].append(a.gapN50())
         GAP_STATS['Total Gap Length'].append(a.gap_length())
       

   # private printer function for what dict is given
   def __print_dict_values(self,dict,list):
      """Prints the dictionary values"""
      for k in list:
         print k,
         for s in dict[k]:
            print DELIMITER, s,
         print

   # public function to print contig stats         
   def print_contig_stats(self):
      """Prints the contig stats"""
      self.__print_dict_values(CON_STATS,C_STATS)

   # public function to print scaffold stats
   def print_scaffold_stats(self):
      """Prints the scaffold stats"""
      self.__print_dict_values(SCAFF_STATS,S_STATS)

   # public function to print gap stats
   def print_gap_stats(self):
      """Prints the gap stats"""
      self.__print_dict_values(GAP_STATS,G_STATS)

   # public function to print info about assembly
   def print_assembly_names(self,assembly):
      """Prints the assembly name, if given"""
      print "Name",
      for a in assembly:
         print DELIMITER, a.get_assembly_name(),
      print

   # public function to print info about assembler
   def print_assembler(self,assembly):
      """Prints the assembler, if given"""
      print "Assembler",
      for a in assembly:
         print DELIMITER, a.get_assembler(),
      print

   # public function to print all assembly stats    
   def print_stats(self):
      """Prints the full assembly stats"""
      self.print_assembly_names(self.assembly)
      self.print_assembler(self.assembly)
      self.print_contig_stats()
      self.print_scaffold_stats()
      self.print_gap_stats()

   def __get_assembly(self):
      return self.assembly
   
   def __get_assembly_names(self):
      names = []
      names.append("Name")
      for a in self.__get_assembly():
         names.append(a.get_assembly_name())
      return names

   def __get_assembler(self):
      assembler = []
      assembler.append("Assembler")
      for a in self.__get_assembly():
         assembler.append(a.get_assembler())
      return assembler

   # private printer function for what dict is given
   def __get_dict_values(self,dict,list):
      """Gets the dictionary values"""
      values = []
      for k in list:
         tmp = []
         tmp.append(k)
         for s in dict[k]:
            tmp.append(s)
         values.append(tmp)
      return values

   def __get_gap_stats(self):
      return self.__get_dict_values(GAP_STATS,G_STATS)

   def __get_scaffold_stats(self):
      return self.__get_dict_values(SCAFF_STATS,S_STATS)

   def __get_contig_stats(self):
      return self.__get_dict_values(CON_STATS,C_STATS)

   def get_assembly_names(self):
      return self.__get_assembly_names()
   
   def get_stats(self):
      data = []
      data.append(self.__get_assembler())
      data += self.__get_contig_stats()
      data += self.__get_scaffold_stats()
      data += self.__get_gap_stats()
      return data
