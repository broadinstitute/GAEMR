'''
Created on Apr 30, 2012

@author: bruce
'''

# Copyright 2012 The Broad Institute, Inc.
# SOFTWARE COPYRIGHT NOTICE
# This software and its documentation are the copyright of the
# Broad Institute, Inc. All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. The Broad Institute is not responsible for its
# use, misuse, or functionality.

import math
#from Bio import SeqRecord

class Kmer:
    """
    Kmer-generating class.
    """
    # 2-bit bases, rc is base ^ 3
    BASES = "ACGT"
    RC = {'A': 'T', 
          'C': 'G',
          'G': 'C',
          'T': 'A'}
        
    def __init__(self, K, rc=True, slide=True):
        """
        Class initializer: 
            K = kmer size
            rc = whether to consider kmers equivalent to their rc
            slide = whether to using sliding-window kmers or step by K
        """
        self.K = K
        self.rc = rc
        self.slide = slide
        


    def kmerize(self, seq):
        """
        Generator to break seq into kmers using 2-bit-per-base integers.
        Somewhat faster (~10%) with k <= 29.
        """

        k = self.K
        n = 0
        kmer = 0
        rckmer = 0
        mask = (1 << (2 * k)) - 1
        shift = 2 * (k - 1)
        for b in seq:
            if b == 'A': v = 0
            elif b == 'C': v = 1
            elif b == 'G': v = 2
            elif b == 'T': v = 3
            else:
                n = kmer = rckmer = 0
                continue
            kmer = ((kmer << 2) & mask) | v
            rckmer = (rckmer >> 2) | ((v ^ 3) << shift)
            n += 1
            if n >= k:
                if self.rc:
                    yield min(kmer, rckmer)
                else:
                    yield kmer
                if not self.slide:
                    n = kmer = rckmer = 0

    def kmer_rc(self, kmer):
        """
        Retrun kmer corresponding to rc of given kmer.
        """

        rckmer = 0
        for i in xrange(self.K):
            # get next 2-bit base
            b = kmer & 3
            # compute its rc
            b_rc = b ^ 3
            # push onto rc kmer
            rckmer = (rckmer << 2) | b_rc
            # shift source kmer to next base
            kmer >>= 2
        return rckmer

                
    def kmerize_string(self, seq):
        """
        Generator to break seq into kmers using strings.
        """

        n = 0
        k = self.K
        empty = 'x' * k
        kmer = empty
        rckmer = empty
        for b in seq:
            r = self.RC.get(b)
            if not r:
                n = 0
                kmer = empty
                rckmer = empty
                continue
            kmer = kmer[1:] + b
            rckmer = r + rckmer[:-1]
            n += 1
            if n >= self.K:
                if self.rc:
                    yield min(kmer, rckmer)
                else:
                    yield kmer
                if not self.slide:
                    n = kmer = rckmer = 0
                
    def kmer_list(self, seq):
        return [k for k in self.kmerize(seq)]
    
    def kmer_to_seq(self, kmer, rc = False):
        seq = ''
        if not rc:
            shift = 2 * (self.K - 1)
            for i in xrange(self.K):
                b = ((kmer >> shift) & 3)
                seq += self.BASES[b]
                kmer <<= 2
        else:
            for i in xrange(self.K):
                b = (kmer & 3) ^ 3
                seq += self.BASES[b]
                kmer >>= 2
        return seq
    
    def kmer_count_seq(self, seq, counts = None):
        if counts is None: counts = {}
        for k in self.kmerize(seq):
            counts[k] = counts.get(k, 0) + 1
        return counts
    
    def kmer_count_seqs(self, seqs):
        counts = {}
        for s in seqs:
            self.kmer_count_seq(s, counts)
        return counts

    def kmer_counts_to_mean(self, counts):
        distinct = 0
        total = 0
        for (kmer, count) in counts.iteritems():
            distinct += 1
            total += count
        if distinct == 0: return 0
        return float(total) / float(distinct)

    def kmer_entropy(self, seqs):
        """
        Computes Shannon entropy of seqs using our K value. Return
        value is entropy bits per base.
        """
        possible = 4.0 ** self.K
        counts = self.kmer_count_seqs(seqs)
        total_kmers = sum(counts.values())
        total = 0.0
        for c in counts:
            p = float(counts[c]) / float(total_kmers)
            total += -p * math.log(p, possible)
        return total



    
    
