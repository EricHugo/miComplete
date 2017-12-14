# Copyright (c) Eric Hugoson.
# See LICENSE for details.

from __future__ import print_function, division
from collections import defaultdict
from itertools import chain
import re

class linkageAnalysis():
    def __init__(self, seqObject, baseName, seqType, proteome, seqstats, 
            hmmMatches, debug=False, q=None):
        self.baseName = baseName
        self.seqObject = seqObject
        self.q = q
        self.seqType = seqType
        fastats, self.seqLength, allLengths, GC = seqstats
        self.proteome = proteome
        self.hmmMatches = hmmMatches
        self.debug = debug
        if seqType == "faa":
            raise TypeError('Sequences for linkage analysis needs to be fna or \
                    gbk')
        with open(self.proteome) as protFile:
            self.pHeaders = set(header for header in protFile 
                    if re.search("^>", header))
            #print(self.pHeaders)

    def get_locations(self):
        """Get locations or start and stop for each matched produced by
        the completeness check, including duplicates. Also hooks in downstream
        functions that eventually calculate the relative linkage of all markers"""
        # prodigal writes ID=X_XX... in proteome, also START # STOP of sequence
        # the ID can be found in .tblout, therefore if match ID => get location
        # Force same format out of .gbk
        #
        #print(self.hmmMatches)
        #
        ## for each matching gene get weight from self.pHeaders and put into
        ## new dict in format dict[hmm] = [[START, STOP], [START, STOP], [...]]
        self.hmmLocations = defaultdict(list)
        for hmm, genes in self.hmmMatches.items():
            for gene in genes:
                for loc in self.pHeaders:
                    if re.search(re.escape(gene[0])+"\s", loc):
                        # convert to int and append to dict[hmm]
                        self.hmmLocations[hmm].append(list(map(int,
                            loc.split('#')[1:3])))
            #print(self.hmmLocations[hmm])
        return self.hmmLocations 

    def find_neighbour_distance(self):
        """For each location of each matched marker the location of start and stop
        are compared to all other locations, lowest value is stored"""
        # for each hmm compare loc(s) to all other locs to find lowest negative
        # value, these are closest neighbours up- and downstream
        try:
            self.hmmLocations
        except AttributeError:
            self.get_locations()
        if self.debug:
            #self.hmmLocations["test"].append([110, 120])
            #self.hmmLocations["test2"].append([510, 720])
            pass
        self.locs = defaultdict(list)
        for hmm, locs in self.hmmLocations.items():
            minFLoc = []
            minRLoc = []
            #print(hmm)
            for loc in locs:
                # nested list comprehension
                # reads locs and compares end of current read to start of all
                # if negative -> adds sequence length to simulate circularity
                forwardL = [[ int(each[0] - loc[1]) if int(each[0] - loc[1] > 0)
                    else int(each[0] - loc[1] + self.seqLength) for each in forw ]
                        for key, forw in self.hmmLocations.items() if not
                        key == hmm ]
                # flatten list
                forwardLFlat = list(chain.from_iterable(forwardL))
                minFLoc.append(min(forwardLFlat))
                reverseL = [[ int(loc[0] - each[1]) if int(loc[0] - each[1] > 0)
                    else int( loc[0] - each[1] + self.seqLength) for each in rev ]
                        for key, rev in self.hmmLocations.items() if not
                        key == hmm ]
                reverseLFlat = list(chain.from_iterable(reverseL))
                minRLoc.append(min(reverseLFlat))
            self.locs[hmm].append(min(minFLoc))
            self.locs[hmm].append(min(minRLoc))
            #print(self.locs[hmm])
        return self.locs

    def calculate_linkage_scores(self):
        """From dict of with smallest distance locations up- and downstream
        the average is calculated, totaled and a relative value is calculated"""
        try:
            self.locs
        except AttributeError:
            self.find_neighbour_distance()
        self.linkageAbsVals = {hmm: (loc[0] + loc[1]) / 2 for (hmm, loc) in
                self.locs.items()}
        self.totalDistance = sum([ linkVal for hmm, linkVal in
            self.linkageAbsVals.items()])
        #print(self.linkageAbsVals)
        #print(self.totalDistance)
        self.linkageRelVals = {hmm: [(linkVal / self.totalDistance)]
                for hmm, linkVal in self.linkageAbsVals.items()}
        #print(self.linkageRelVals)
        # Send results to function with lock to write to single file
        # Once all results are there average and make boxplot
        # Write resulting weights to .weights file
        return self.linkageRelVals

