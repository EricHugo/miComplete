# Copyright (c) Eric Hugoson.
# See LICENSE for details.

from __future__ import print_function, division
from distutils import spawn
from collections import defaultdict
from termcolor import cprint
import sys
import math
import subprocess
import re
import os


class calcCompleteness():
    def __init__(self, fasta, baseName, argv, linkage=False):
        """Initializes the basic variables of the completeness search"""
        self.baseName = baseName
        self.evalue = "-E %s" % (str(argv.evalue))
        self.tblout = "%s.tblout" % (self.baseName)
        self.hmms = argv.hmms
        self.fasta = fasta
        self.linkage = linkage
        self.argv = argv
        print("Starting completeness for " + fasta, file=sys.stderr)
        self.hmmNames = set({})
        with open(self.hmms) as hmmfile:
            for line in hmmfile:
                if re.search('^NAME', line):
                    name = line.split(' ')
                    self.hmmNames.add(name[2].strip())
        if self.argv.debug:
            print(self.hmmNames)

    def hmm_search(self):
        """Runs hmmsearch using the supplied .hmm file and produces a table
        out, also hooks in get_completeness() and ultimately returns found
        markers, duplicated markers and total number of markers"""
        hmmsearch = ["hmmsearch", self.evalue, "--tblout", self.tblout,
                self.hmms, self.fasta]
        if self.argv.debug:
            print(hmmsearch)
        if sys.version_info > (3, 4):
            compProc = subprocess.run(hmmsearch, stdout=subprocess.DEVNULL)
            errcode = compProc.returncode
        else:
            errcode = subprocess.call(hmmsearch, stdout=open(os.devnull, 'wb'),
                    stderr=subprocess.STDOUT)
        if errcode > 0:
            cprint("Warning:", 'red', end=' ')
            print("Error thrown by HMMER, is %s empty?" % self.fasta)
            return 0, 0, 0
        if self.linkage:
            self.get_completeness()
            return self.filledHmms
        self.get_completeness()
        return self.numHmms, self.numDupHmms, len(self.hmmNames)

    def get_completeness(self):
        """Reads the out table of hmmer to find which hmms are present, and
        which are duplicated"""
        self.hmmMatches = defaultdict(list)
        self.seenHmms = set()
        # gather gene name and evalue in dict by key[hmm]
        for hmm in self.hmmNames:
            with open(self.tblout) as hmmTable:
                for foundHmm in hmmTable:
                    if re.match("#$", foundHmm):
                        break
                    if re.search("^" + hmm + "$", foundHmm.split()[2]):
                        #slice notation gathers column 0 and 4
                        self.hmmMatches[hmm].append(foundHmm.split()[0:5:4])
                        self.seenHmms.add(hmm)
        if self.argv.debug:
            print(self.hmmMatches)
            print(len(self.hmmMatches))
        self.filledHmms = defaultdict(list)
        # section can be expanded to check for unique gene matches
        for hmm, geneMatches in self.hmmMatches.items():
            # sort by lowest eval to fill lowest first
            for gene in sorted(geneMatches, key=lambda ev: float(ev[1])):
                if hmm not in self.filledHmms:
                    self.filledHmms[hmm].append(gene)
                elif float(gene[1]) < pow(float(self.filledHmms[hmm][0][1]), 1/2):
                    self.filledHmms[hmm].append(gene)
        self.numHmms = len(self.filledHmms)
        self.dupHmms = [hmm for hmm, genes in self.filledHmms.items()
                if len(genes) > 1]
        self.numDupHmms = len(self.dupHmms)
        if self.argv.hlist and not self.linkage:
            self.print_hmm_lists()
        return

    def print_hmm_lists(self):
        """Prints the contents of found, duplicate and and not found markers"""
        hlistName = "%s_hmms.list" % (self.baseName)
        with open(hlistName, 'w+') as seenList:
            for eachHmm in self.seenHmms:
                seenList.write("%s\n" % eachHmm)
        dupListName = "%s_hmms_duplicate.list" % (self.baseName)
        with open(dupListName, 'w+') as dupList:
            for eachDup in self.dupHmms:
                dupList.write("%s\n" % eachDup)
        missingListName = "%s_hmms_missing.list" % (self.baseName)
        with open(missingListName, 'w+') as missingList:
            for hmm in self.hmmNames:
                if hmm not in self.seenHmms:
                    missingList.write("%s\n" % hmm)

    def attribute_weights(self, numHmms):
        """Using the markers found and duplicates from get_completeness(), and
        provided weights, adds up weight of present and duplicate markers"""
        weightedComplete = 0
        weightedRedun = 0
        for hmm in self.seenHmms:
            with open(self.argv.weights, 'r') as weights:
                for eachWeight in weights:
                    if re.match(hmm + "\s", eachWeight):
                        weightedComplete += float(eachWeight.split()[1])
        for hmm in self.dupHmms:
            with open(self.argv.weights, 'r') as weights:
                for eachWeights in weights:
                    if re.match(hmm + "\t", eachWeights):
                        weightedRedun += float(eachWeights.split()[1])
        weightedRedun = round(((weightedRedun + weightedComplete) /
            weightedComplete), 3)
        weightedComplete = round(weightedComplete, 3)
        return weightedComplete, weightedRedun


