
# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re
import os
from RunCommand import RunCommand
import PlatformConstant as pc

constant = pc.PlatformConstant()

class Taxonomy:
    """This class helps with BLAST nodes and names taxonomy files"""

    def __init__(self, nodes_db=None, names_db=None, blast_db=None):
        self.nodes_obj = Nodes(nodes_db)
        self.names_obj = Names(names_db)
        self.blastdbcmd_obj = BlastDbCmd(blast_db)

    def __get_nodes(self):
        return self.nodes_obj.get_nodes_data()

    def __get_names(self):
        return self.names_obj.get_names_data()

    def have_nodes(self):
        return self.__get_nodes() != 'None'

    def have_names(self):
        return self.__get_names() != 'None'

    # private function to get the lineage for this taxonomic id
    def get_taxonomy_string(self, id, is_gi_number=False):
        """Returns lineage for a particular tax id."""
        tax_id = id
        if is_gi_number:
            tax_id = self.get_tax_id_from_gi(id)
        tmp = []
        first = 1
        nodes = self.__get_nodes()
        names = self.__get_names()

        if nodes and names:
            while nodes[tax_id]["parent"] and nodes[tax_id]["parent"] != tax_id and names[tax_id] != "root":
                if nodes[tax_id]["rank"] == 'no rank':
                    if first:
                        tmp.append("common_name=" + names[tax_id])
                    else:
                        tmp.append("domain=" + names[tax_id])
                else:
                    tmp.append(nodes[tax_id]["rank"] + "=" + names[tax_id])
                    first = 0
                tax_id = nodes[tax_id]["parent"]

            if tmp:
                tmp.reverse()
                return ''.join([i + ";" for i in tmp]).rstrip(';')

        return "domain=NOT_FOUND"

    # private function to just get genus
    def get_taxonomic_level(self, id, level, is_gi_number=False):
        """Returns genus name for a given tax id."""
        nodes = self.__get_nodes()
        names = self.__get_names()
        if nodes and names:

            tax_id = id
            if is_gi_number:
                tax_id = self.get_tax_id_from_gi(gi)

            if tax_id not in nodes or tax_id not in names:
                return "None"
            while nodes[tax_id]["parent"] and nodes[tax_id]["parent"] != tax_id and names[tax_id] != "root":
                if nodes[tax_id]["rank"] == level:
                    break
                tax_id = nodes[tax_id]["parent"]

            return names[tax_id]
        return None

    def get_domain(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'domain', is_gi_number)

    def get_kingdom(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'kingdom', is_gi_number)

    def get_phylum(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'phylum', is_gi_number)

    def get_class(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'class', is_gi_number)

    def get_order(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'order', is_gi_number)

    def get_family(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'family', is_gi_number)

    def get_genus(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'genus', is_gi_number)

    def get_species(self, id, is_gi_number=False):
        return self.get_taxonomic_level(id, 'species', is_gi_number)

    def get_tax_id_from_gi(self, gi):
        return self.blastdbcmd_obj.get_tax_id_from_gi(gi)

    def get_gi_tax_lookup(self, gi_list):
        return self.blastdbcmd_obj.get_gi_tax_lookup(gi_list)

    def get_gi_tax_name_lookup_by_level(self, gi_tax_lookup, level='genus'):
        gi_tax_name_lookup = {}
        for gi in gi_tax_lookup:
            gi_tax_name_lookup[gi] = self.get_taxonomic_level(gi_tax_lookup[gi], level)
        return gi_tax_name_lookup
    
class Nodes:

    def __init__(self, nodes_db=None):
        self.nodes = None
        if nodes_db:
            self.nodes = self.__load_nodes(nodes_db)

    # private function to read the blast db nodes taxdump file
    def __load_nodes(self, db):
        """Loads nodes dict with parent and rank values."""
        nodes = {}
        try:
            f=open(db)
            for lines in f.readlines():
                nodes_data = lines.rstrip('\n').split('\t|\t')
                key = nodes_data[0]
                nodes[key] = {}
                nodes[key]["parent"] = nodes_data[1]
                nodes[key]["rank"] = nodes_data[2]
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1
        return nodes

    def get_nodes_data(self):
        return self.nodes

class Names:

    def __init__(self, names_db):
        self.names = None
        if names_db:
            self.names = self.__load_names(names_db)

    # private function to read the blast db names taxdump file
    def __load_names(self, db):
        """Loads names dict with tax id and names."""
        names = {}
        try:
            f=open(db)
            for lines in f.readlines():
                if re.search("scientific",lines):
                    names_data = lines.rstrip('\n').split('\t|\t')
                    re.sub("=","-",names_data[1])
                    names[names_data[0]] = names_data[1]
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1
        return names

    def get_names_data(self):
        return self.names


BLAST_DB_OUTPUT = ['gi','accession','ordinal id','sequence title','length','tax id','common name','scientific name']

class BlastDbCmd:

    def __init__(self, blast_db=None, gi_list=None):
        self.blast_db = blast_db
        self.blast_hit_info = {}
        if gi_list:
            self.blast_hit_info = self.__check_gi_list(gi_list)

    def __get_blastdb(self):
        return self.blast_db

    def __get_gi_string(self, gi_nums):
        """Returns string of gi numbers."""
        # Deprecated 'cuz we had cases of too many gi nums
        gi_string = ''.join([gi_nums[i] + "," for i in range(int(len(gi_nums)) - 1)])
        gi_string += gi_nums[-1]
        return gi_string

    # private function to get the command to run
    def __build_command(self, gi_file):
        """Returns a blastdbcmd command string."""
        outfmt_string = ''.join(['%' + i + constant.TABLE_DELIMITER for i in ['g','a','o','t','l','T','L','S']])
        return[constant.BLASTDBCMD,"-db",self.__get_blastdb(),"-outfmt",outfmt_string,"-entry_batch",gi_file]

    def __populate_blast_hit_info(self, data):
        for i in data:
            if i: # just in case there is a blank line
                tax_data = i.split(constant.TABLE_DELIMITER)
                if tax_data[0] not in self.blast_hit_info:
                    self.blast_hit_info[tax_data[0]] = {}
                for j in xrange(1,len(tax_data) - 1):
                    self.blast_hit_info[tax_data[0]][BLAST_DB_OUTPUT[j]] = tax_data[j]

    def __check_gi_list(self, gi_list):
        tmp_gi_list = []
        for gi in gi_list:
            if not self.__have_gi_info(gi):
                tmp_gi_list.append(gi)
        if tmp_gi_list:
            try:
                tmp_gi_file = 'tmp.gi_list'
                out = open(tmp_gi_file, 'w')
            except:
                print "Can't open tmp gi list file, " + tmp_gi_file
                sys.exit(-1)
            for i in tmp_gi_list:
                out.write(i + '\n')
            out.close()
            command = self.__build_command(tmp_gi_file)
            if command:
                rc = RunCommand(command)
                cmd_data = rc.run_command(1,0).split('\n')
                self.__populate_blast_hit_info(cmd_data)
            else:
                print "Unable to create blastdbcmd.  Exiting."
                sys.exit(-1)
            os.remove(tmp_gi_file)

    def __have_gi_info(self, gi):
        return gi in self.blast_hit_info

    def __get_info_helper(self, gi, key=None):
        if not self.__have_gi_info(gi):
            self.__check_gi_list([gi])
        return self.blast_hit_info[gi][key]

    def get_accession_from_gi(self, gi):
        return self.__get_info_helper(gi, 'accession')

    def get_ordinal_id_from_gi(self, gi):
        return self.__get_info_helper(gi, 'ordinal id')
    
    def get_sequence_title_from_gi(self, gi):
        return self.__get_info_helper(gi, 'sequence title')
        
    def get_length_from_gi(self, gi):
        return self.__get_info_helper(gi, 'length')

    def get_tax_id_from_gi(self, gi):
        return self.__get_info_helper(gi, 'tax id')

    def get_common_name_from_gi(self, gi):
        return self.__get_info_helper(gi, 'common name')

    def get_scientific_name_from_gi(self, gi):
        return self.__get_info_helper(gi, 'scientific name')
    
    def get_gi_tax_lookup(self, gi_list=None):
        """Loads dict with gi to tax id values."""
        if gi_list:
            self.__check_gi_list(gi_list)
        gi_tax_dict = {}
        for i in self.blast_hit_info.keys():
            gi_tax_dict[i] = self.blast_hit_info[i]['tax id']
        return gi_tax_dict

    # Output format, where the available format specifiers are:
   # %f means sequence in FASTA format
   # %s means sequence data (without defline)
   # %a means accession
   # %g means gi
   # %o means ordinal id (OID)
   # %t means sequence title
   # %l means sequence length
   # %T means taxid
   # %L means common taxonomic name
   # %S means scientific name
   # %P means PIG
