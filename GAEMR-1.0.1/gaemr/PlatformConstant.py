#!/usr/bin/env python

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import os

class PlatformConstant(object):
    """A unified class to contain all platform constants"""

    MIN_GAP_SIZE=10
    MIN_CONTIG_SIZE=200
    MIN_SCAFFOLD_SIZE=1000
    QUAL_MAX=93
    MAX_INSERT_SIZE=50000
    TABLE_DELIMITER=" | "
    MUMMER_PATH="/broad/software/groups/gtba/software/mummer_3.23-64bit/" # "/path/to/mummer/package/" # http://mummer.sourceforge.net/
    BLAST_DIR="/broad/software/groups/gtba/software/ncbi-blast-2.2.25+/bin/" # "/path/to/blast+/package/" # http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download


    SAMTOOLS="/broad/software/groups/gtba/software/samtools_0.1.18/bin/samtools" # "/path/to/samtools/samtools" # http://samtools.sourceforge.net/
    PICARD="/seq/software/picard/current/bin/" # "/path/to/public/picard/package/" # http://sourceforge.net/projects/picard/
    PICARD_TMP_DIR="."
    BLAST_NT="/broad/data/blastdb/nt/nt" # "nt"
    BLAST_UNIVEC="/gsap/assembly_analysis/databases/UniVec/UniVec" # "/path/to/UniVec/db/db" # ftp://ftp.ncbi.nih.gov/pub/UniVec/
    BLAST_rRNA="/gsap/assembly_analysis/databases/NCBI_rRNA/ncbi_rRNA" # "/path/to/rRNA/db/db" # Manually curated
    BLAST_MITOGCONTAM="/gsap/assembly_analysis/databases/mitogcontam/mitogcontam" # "/path/to/mitogcontam/db/db" # ftp://ftp.ncbi.nih.gov/pub/kitts/gcontam1.gz & ftp://ftp.ncbi.nih.gov/refseq/release/mitochondrion/*.f*a.gz
    BLAST_NODES="/broad/data/taxonomy/taxdump/nodes.dmp" # "/path/to/taxdump/nodes.dmp" # ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    BLAST_NAMES="/broad/data/taxonomy/taxdump/names.dmp" # "/path/to/taxdump/names.dmp" # ftp://ftp.ncbi.nih.gov/pub/taxonomy/
    RNAMMER="/seq/annotation/bio_tools/rnammer/current/rnammer" # "/path/to/rnammer/package/rnammer" # http://www.cbs.dtu.dk/cgi-bin/sw_request?rnammer
    RDP="/broad/software/groups/gtba/software/rdp_classifier_2.4/rdp_classifier-2.4.jar" # "/path/to/rdp_classifier_2.4/package/rdp_classifier-2.4.jar" # http://sourceforge.net/projects/rdp-classifier/

    BLASTDBCMD=BLAST_DIR + "blastdbcmd"
    MAKEBLASTDB=BLAST_DIR + "makeblastdb"
    NUCMER=MUMMER_PATH + "nucmer"
    PROMER=MUMMER_PATH + "promer"
    SHOWTILING=MUMMER_PATH + "show-tiling"
    MUMMERPLOT=MUMMER_PATH + "mummerplot"

    BWA="/seq/software/picard/current/3rd_party/bwa" # "/path/to/bwa/package/" # http://bio-bwa.sourceforge.net/
    GAEMR="/gsap/assembly_analysis/GAEMR/bin/" # /your/local/install/GAEMR
    BOWTIE="/broad/software/free/Linux/redhat_5_x86_64/pkgs/bowtie2_2.0.0-beta5/" # "/path/to/bowtie/package/" # http://bowtie-bio.sourceforge.net/index.shtml
    PLOT_COLORS=['green','blue','red','orange','magenta','black','grey']

    #detect user, environment and populate field
    def __init__(self):
        """populate class fields"""
        self.user_name = os.getenv("USER")
        self.operating_system = os.getenv("OSTYPE")
        
