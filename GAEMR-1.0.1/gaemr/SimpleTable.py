#
# Class for storing and formatting table data
#

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import sys
import re
import PlatformConstant as pc
constant = pc.PlatformConstant()

class SimpleTable():
    """ This class represents a manipulation of Fastq files. """

    def __init__(self, headers, data, title, empty_data_default="NA"):
        self.title = title
        self.headers,self.data = self.__check_data_lengths(headers,data,empty_data_default)
        # need to check length of headers and data
        # fill with default if < len(headers); die if > len(headers)

    def __check_data_lengths(self, headers, data, default):
        header_len = len(headers)

        for i in xrange(len(data)):
            data_len = len(data[i])
            assert data_len <= header_len
            for j in range(len(data[i]),header_len):
                data[i].append(default)
        return headers, data

    def add_row(self,row_list):
        data = self.get_data()
        data.append(row_list)
        
    def get_headers(self):
        return self.headers

    def get_data(self):
        return self.data

    def get_title(self):
        return self.title

    def __add_commas_to_nums(self, num):
        formatted_string = ''
        #print "NUM",\
        #print num
        data = num.split(' ')
        #print "DATA:",\
        #print data
        for i in xrange(len(data)):
            value = data[i]
            if re.match('\d+', str(value)) and not re.search('\D+', str(value)):
                if re.search('\.', str(value)):
                    value = "{:,.2f}".format(float(value))
                else:
                    value = "{:,}".format(int(value))
            formatted_string += value + ' '
        return formatted_string[:-1]

    #delete this comment in a few days after no dependency problems are confirmed (9/27)
    #def to_table(self,delimiter=constant.TABLE_DELIMITER):
    def to_table(self,delimiter=' | '):
        title = "#TITLE: " + str(self.get_title()) + '\n'
        header_line = "#"
        header_line += ''.join([str(i) + delimiter for i in self.get_headers()])
        header_line = header_line.rstrip(delimiter) + '\n'
        data = self.get_data()
        data_line = ''
        for i in xrange(len(data)):
            data_line += ''.join([self.__add_commas_to_nums(str(data[i][j])) + delimiter for j in xrange(len(data[i]))])
            data_line = data_line.rstrip(delimiter) + '\n'
        return title + header_line + data_line

    def __get_html_title(self):
        return ' <h1>' + self.get_title() + '</h1>'

    def __get_html_header(self):
        headers = self.get_headers()
        html_header_string = '   <tr>\n'
        html_header_string += ''.join(['    <th width=\"25%\">' + str(i) + '</th>\n' for i in headers])
        html_header_string += '   </tr>\n'
        return html_header_string

    def __get_html_data(self):
        data = self.get_data()
        html_string = ''
        for i in xrange(len(data)):
            for j in xrange(len(data[i])):
                if j == 0:
                    html_string += '   <tr>\n    <th>' + self.__add_commas_to_nums(str(data[i][j])) + '</th>\n'
                    continue
                html_string += '    <td>' + self.__add_commas_to_nums(str(data[i][j])) + '</td>\n'
            html_string += '   </tr>\n'
        return html_string
        
    def to_html(self):
        #html_string = '<div name=\"' + self.__get_html_title() + '\">\n'
        #html_string += self.__get_html_title() + '\n'
        html_string = '  <table id=\"' + self.get_title() + '\">\n'
        #html_string += '  <table>\n'
        html_string += self.__get_html_header()
        html_string += self.__get_html_data()
        html_string += '  </table>\n'
        #html_string += '</div>\n'
        return html_string

    def print_output(self,header=None,html=True):
        text_table = sys.stdout
        html_table = sys.stdout
        if header:
            try:
                text_table = open(header + '.table.txt','w')
                if html:
                    html_table = open(header + '.table.html','w')
            except IOError as (errno,strerror):
                print "I/O Error({0}): {1}".format(errno,strerror)
                return -1
        text_table.write(self.to_table())
        if html:
            html_table.write(self.to_html())
