# Copyright (c) Eric Hugoson.
# See LICENSE for details.


from __future__ import print_function, division
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

class parseSeqStats():
    def __init__(self, seq, baseName, seqType):
        """Initialize generators for headers as well as sequence"""
        if seqType == "fna" or seqType == "faa":
            self.seqType = "fasta"
        elif re.match("(gb.?.?)|genbank", seqType):
            self.seqType = "genbank"
        else:
            self.seqType = seqType
        self.seqHeaders = [seqi.id for seqi in SeqIO.parse(seq, self.seqType)]
        self.seqFasta = [seqi.seq for seqi in SeqIO.parse(seq, self.seqType)]
        self.baseName = baseName
        self.sequence = seq

    def get_stats(self, totalLength, allLengths):
        """Finds assembly stats if applicable"""
        length = []
        N50, N90 = "", ""
        for contig in sorted(allLengths, reverse=True):
            length.append(contig)
            if sum(length) >= totalLength * 0.9:
                N90 = contig
                L90 = len(length)
                if N50:
                    break
            if sum(length) >= totalLength / 2:
                if N50:
                    continue
                N50 = contig
                L50 = len(length)
                if N90:
                    break
        return N50, L50, N90, L90

    def get_length(self):
        allLengths, totalFasta = [], []
        for fasta in self.seqFasta:
            allLengths.append(len(fasta))
            totalFasta.append(str(fasta))
        if self.seqType == "faa":
            GCcontent = 0.0
        else:
            # implement memory-efficient method here
            GCcontent = round(GC(''.join(str(totalFasta))), 2)
        seqLength = sum(allLengths)
        return seqLength, allLengths, GCcontent

