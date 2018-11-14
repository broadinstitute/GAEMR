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
import re
from Bio.Blast.Applications import *
from optparse import OptionParser
from gaemr.RunCommand import RunCommand
import gaemr.PlatformConstant as pc

constant = pc.PlatformConstant()

parser = OptionParser(usage="usage: %prog [options] <db file> <query file>")

parser.add_option('--evalue', '-e', action="store", default="1e-5", type='string', dest="evalue",
    help='E-value cutoff (default=%default)')
parser.add_option('--outformat', '-f', action="store", default="5", type='int', dest="outfmt",
    help='Output format for blast results (default=%default)')
parser.add_option('--output', '-o', action="store", default=None, type='string', dest="output",
    help='Blast output file name (default=%default)')
parser.add_option('--threads', '-t', action="store", default="4", type='int', dest="threads",
    help='Number of threads to use for blastn (default=%default)')
parser.add_option('--max_targets', '-m', action="store", default="10", type='int', dest="max_targets",
    help='The maximum number of target hits to report (default=%default)')
parser.add_option('--task', '-b', action="store", default="blastn", type='string', dest="task",
    help='The blast task to run [blastn,blastn-short,megablast,\
                  dc-megablast,blastp,blastx,tblastn,tblastx] (default=%default)')
parser.add_option('--vecscreen', '-v', action="store_true", default=False, dest="vec_screen",
    help='Run blastn mimicing VecScreen parameters (default=%default)')
parser.add_option('--rRNA', '-r', action="store_true", default=False, dest="rRNA_screen",
    help='Run megablast mimicing NCBI rRNA screen parameters (default=%default)')
parser.add_option('--ncbi', '-n', action="store_true", default=False, dest="ncbi_screen",
    help='Run megablast mimicing NCBI gcontam parameters (default=%default)')
parser.add_option('--makeblastdb', '-d', action="store_true", default=False, dest="make_blastdb",
    help='Format the blast db using makeblastdb (default=%default)')
parser.add_option('--is_protein', '-p', action="store_true", default=False, dest="is_protein",
    help='BLAST db to make is a protein db (default=%default)')

(options, args) = parser.parse_args()

TASKS = {'blastn': NcbiblastnCommandline,
         'blastn-short': NcbiblastnCommandline,
         'megablast': NcbiblastnCommandline,
         'dc-megablast': NcbiblastnCommandline,
         'blastp': NcbiblastpCommandline,
         'blastx': NcbiblastxCommandline,
         'tblastn': NcbitblastnCommandline,
         'tblastx': NcbitblastxCommandline
}

#- NcbiblastpCommandline - Protein-Protein BLAST
#- NcbiblastnCommandline - Nucleotide-Nucleotide BLAST
#- NcbiblastxCommandline - Translated Query-Protein Subject BLAST
#- NcbitblastnCommandline - Protein Query-Translated Subject BLAST
#- NcbitblastxCommandline - Translated Query-Protein Subject BLAST


# check inputs or die
if len(args) != 2:
    parser.error("Must supply db and query fasta files.")
    sys.exit(-1)

if not options.output:
    parser.error("Must supply output file with --output.")
    sys.exit(-1)

### VECSCREEN BLAST+ BLASTN OPTIONS:
# reward = 1
# penalty = -5
# gapopen = 3
# gapextend = 3
# searchsp = 1.75e12 -- needs to be int value
# 1750000000000
# dust = yes


###NCBI gcontam Screen options:
# task = megablast
# dust = yes
# soft_masking = true
# perc_identity = 90
# lcase_masking = ''


###NCBI rRNA Screen options
# blastall options:
## megablast -A 120 -D 3 -E 2 -e 1e-9 -F "m;D" -G 4 -g T -M 5000000 -n T -q -4 -R -r 3 -t 18 -U T -W 12 -X 20 -p 95
# task = megablast
# window_size = 120
# db_gencode = 3
# gapextend = 2
# evalue = 1e-9
# dust = yes
# soft_masking = true
# gapopen = 4
# (-G not needed)
# matrix	= 5000000
# no_greedy = ''
# penalty = -4
# in_pssm
# reward = 3
# max_intron_length = 18
# lcase_masking = ''
# word_size = 12
# xdrop_gap = 20
# perc_identity = 95





def main():
    # set up command
    cmd = None
    task = options.task

    if options.make_blastdb:
        type = 'nucl'
        if options.is_protein:
            type = 'prot'
        cmd_list = [constant.MAKEBLASTDB, '-in', args[0], '-dbtype', type]
        rc = RunCommand(cmd_list)
        print "Running command:  " + rc.get_command() + '\n'
        rc.run_command()

    if options.vec_screen:
        cmd = NcbiblastnCommandline(query=args[1], db=args[0], evalue=700, outfmt=options.outfmt, reward=1,
            penalty=-5, gapopen=3, gapextend=3, dust='yes', searchsp=1750000000000,
            out=options.output, task=task, num_threads=options.threads)
    elif options.ncbi_screen:
        cmd = NcbiblastnCommandline(query=args[1], db=args[0], outfmt=options.outfmt, dust='yes',
            perc_identity=90, lcase_masking='', task='megablast',
            out=options.output, num_threads=options.threads)
        #soft_masking='true',
    elif options.rRNA_screen:
        cmd = NcbiblastnCommandline(query=args[1], db=args[0], outfmt=options.outfmt, dust='yes',
            perc_identity=95, lcase_masking='', task='megablast',
            out=options.output, num_threads=options.threads, evalue=1e-9, window_size=120, gapextend=2, gapopen=4,
            no_greedy='', penalty=-4, reward=3, word_size=12, xdrop_gap=20)
            #Options not supported: in_pssm='',soft_masking='true', matrix=5000000, max_intron_length=18, db_gencode=3,
    else:
        if options.task in TASKS:
            program = TASKS[options.task]
            if re.search("Ncbiblastn", str(program)):
                cmd = program(query=args[1], db=args[0], evalue=options.evalue,
                    outfmt=options.outfmt, out=options.output,
                    num_threads=options.threads, max_target_seqs=options.max_targets,
                    task=task)
            else:
                cmd = program(query=args[1], db=args[0], evalue=options.evalue,
                    outfmt=options.outfmt, out=options.output,
                    num_threads=options.threads, max_target_seqs=options.max_targets)
        else:
            print "Unrecognized blast task, " + options.task
            sys.exit(-1)

    print "Running BLAST command:  " + str(cmd) + '\n'
    out, err = cmd()

    return 0

if __name__ == "__main__":
    sys.exit(main())
    
