
# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import re

class RnammerFile:
    """This class helps with rnammer output fasta file format"""
    def __init__(self,file):
        self.file = file

        f_header_list = self.__get_headers_from_file()
        self.molecules = self.__get_rna_data(f_header_list)

    # private function to get info from fasta header lines
    def __get_headers_from_file(self):
        headers = []
        try:
            f = open(self.file,'r')
            for line in f.readlines():
                if re.match(">",line):
                    headers.append(re.sub(">","",line).strip('\n'))
            f.close()
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1
        return headers

    # private function to get the info out of the header
    def __get_rna_data(self,headers):
        molecules = {}

        # go through all headers
        for i in headers:
            # process the components
            data = i.split(' ')
            name_string = data[0]
            molecule = re.sub("/\w+=","",data[1]).rstrip('_rRNA')
            contig_info = data[0].split('_')
            contig = ''.join([str(x) + '_' for x in contig_info[1:-2]])[:-1]
            dir = re.sub("DIR","",contig_info[-1])
            start,end = contig_info[-2].split('-')
            score = re.sub("/\w+=","",data[2])

            # build molecule object
            m = Molecule(molecule, contig, start, end, dir, score, name_string)

            if not molecule in molecules:
                molecules[molecule] = []
            molecules[molecule].append(m)

        return molecules

    def get_molecules(self):
        return self.molecules

class RDPFile:
    """This class helps with rdp output file format"""
    def __init__(self, file, molecules, lineage=None, score=None):
        self.file = file
        self.rdp_info, self.order_list = self.__get_rdp_info()
        self.__update_taxonomy_info(molecules, lineage, score)

    # private function to read in rdp file
    def __get_rdp_info(self):
        lines = []

        try:
            f = open(self.file,'r')
            lines = f.readlines()
            f.close()
        except IOError as (errno,strerror):
            print "I/O Error({0}): {1}".format(errno,strerror)
            return -1

        return self.__populate_rdp_dict(lines)

    # private function to build dict of gene name, lineage group, lineage, score
    def __populate_rdp_dict(self,lines):
        rdp = {}

        order_list = []
        first = 1
        for i in lines:
            data = i.split('\t')
            name = data[0]
            rdp[name] = {}
            for j in range(5,int(len(data)),3):
                if not data[j + 1] in rdp[name]:
                    rdp[name][data[j + 1]] = {}
                rdp[name][data[j + 1]]['lineage'] = data[j]
                rdp[name][data[j + 1]]['score'] = float(data[j + 2])
                if first:
                    order_list.append(data[j + 1])
            first = 0

        return rdp, order_list

    # private function to check score
    def __check_score(self, score, check):
        return score > check

    # private function to build up the taxonomic string for the molecule
    def __build_tax_string(self, rdp, order, lineage, score):
        tax_string = ''
        for i in order:
            # had to split it in two lines b/c of weird insistence by python
            # that score was a tuple, though I casted it as string
            if i in rdp:
                tax_string += rdp[i]['lineage'] + ":  "
                tax_string += str(rdp[i]['score']) + "; "

            # if we are at user specified lineage
            if i == lineage:
                # if we have a score
                if score:
                    # if we pass the score threshold, print out string, or print filtered string
                    if self.__check_score(rdp[i]['score'], score):
                        return re.sub("; $","",tax_string)
                    else:
                        return "Filtered taxonomy: " + i + ", " + rdp[i]['lineage'] + ", has score, " + str(rdp[i]['score'])
                return re.sub("; $","",tax_string)
        return re.sub("; $","",tax_string)

    # private function to include taxonomy in the molecule object
    def __update_taxonomy_info(self,molecules, lineage, score):
        rdp = self.get_rdp_info()
        order_list = self.__get_order_list()

        for i in rdp:
            for mol in molecules:
                for j in range(int(len(molecules[mol]))):
                    m_obj = molecules[mol][j]
                    if m_obj.is_same_name(i):
                        m_obj.add_genus(rdp[i]["genus"]['lineage'])
                        m_obj.add_taxonomy_info(self.__build_tax_string(rdp[i], order_list, lineage, score))
                        break

    # private function to get taxonomic order
    def __get_order_list(self):
        return self.order_list

    def get_rdp_info(self):
        return self.rdp_info

class Molecule:
    """This class represents a ribosomal gene molecule"""
    def __init__(self, molecule, contig, start, end, direction, score, name_string):
        self.molecule = molecule
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.direction = direction
        self.score = score
        self.name_string = name_string
        self.taxonomy = None
        self.genus = None

    # public function to update taxonomy stuff
    def add_taxonomy_info(self,taxonomy_string):
        self.taxonomy = taxonomy_string

    # public function to update genus
    def add_genus(self,genus):
        self.genus = genus

    # public getter functions
    def get_genus(self):
        return self.genus

    def get_taxonomy(self):
        if self.taxonomy:
            return self.taxonomy
        return "Unknown"

    # public function to see if the molecule is a certain type (e.g. 16s)
    def is_molecule_type(self,type):
        return self.molecule == type

    # return type of gene
    def get_molecule(self):
        return self.molecule

    # return length of molecule
    def get_length(self):
        return str(self.end - self.start) + " bp"

    # public getter functions (score is rnammer score)
    def get_score(self):
        return self.score

    def get_contig(self):
        return self.contig

    def get_direction(self):
        return self.direction

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def get_name_string(self):
        return self.name_string

    def is_same_name(self, check_name):
        return self.name_string == check_name
