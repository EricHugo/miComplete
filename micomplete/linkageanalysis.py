# Copyright (c) Eric Hugoson.
# See LICENSE for details.

"""Returns the relative distances between all identified hmms within
genome"""

from __future__ import print_function, division
from collections import defaultdict
from itertools import chain
import re

class linkageAnalysis():
    def __init__(self, seq_object, base_name, seq_type, proteome, seqstats,
                 hmm_matches, debug=False, q=None):
        self.base_name = base_name
        self.seq_object = seq_object
        self.q = q
        self.seq_type = seq_type
        _, self.seq_length, _, _ = seqstats
        self.proteome = proteome
        self.hmm_matches = hmm_matches
        self.debug = debug
        self.hmm_locations = defaultdict(list)
        self.locs = defaultdict(list)
        if seq_type == "faa":
            raise TypeError('Sequences for linkage analysis needs to be fna or \
                    gbk')
        with open(self.proteome) as prot_file:
            self.p_headers = set(header for header in prot_file
                                 if re.search("^>", header))
            #print(self.p_headers)

    def get_locations(self):
        """Get locations or start and stop for each matched produced by
        the completeness check, including duplicates. Also hooks in downstream
        functions that eventually calculate the relative linkage of all markers"""
        # prodigal writes ID=X_XX... in proteome, also START # STOP of sequence
        # the ID can be found in .tblout, therefore if match ID => get location
        # Force same format out of .gbk
        #
        #print(self.hmm_matches)
        #
        ## for each matching gene get weight from self.p_headers and put into
        ## new dict in format dict[hmm] = [[START, STOP], [START, STOP], [...]]
        for hmm, genes in self.hmm_matches.items():
            for gene in genes:
                for loc in self.p_headers:
                    if re.search(re.escape(gene[0])+"\s", loc):
                        # convert to int and append to dict[hmm]
                        self.hmm_locations[hmm].append(list(map(int,
                                                      loc.split('#')[1:3])))
            #print(self.hmm_locations[hmm])
        return self.hmm_locations

    def find_neighbour_distance(self):
        """For each location of each matched marker the location of start and stop
        are compared to all other locations, lowest value is stored"""
        # for each hmm compare loc(s) to all other locs to find lowest negative
        # value, these are closest neighbours up- and downstream
        try:
            self.hmm_locations
        except AttributeError:
            self.get_locations()
        if self.debug:
            #self.hmm_locations["test"].append([110, 120])
            #self.hmm_locations["test2"].append([510, 720])
            pass
        for hmm, locs in self.hmm_locations.items():
            min_floc = []
            min_rloc = []
            #print(hmm)
            for loc in locs:
                # nested list comprehension
                # reads locs and compares end of current read to start of all
                # if negative -> adds sequence length to simulate circularity
                forward_l = [[int(each[0] - loc[1]) if int(each[0] - loc[1] > 0)
                              else int(each[0] - loc[1] + self.seq_length) for
                              each in forw]
                             for key, forw in self.hmm_locations.items() if not
                             key == hmm]
                # flatten list
                forward_l_flat = list(chain.from_iterable(forward_l))
                min_floc.append(min(forward_l_flat))
                reverse_l = [[int(loc[0] - each[1]) if int(loc[0] - each[1] > 0)
                              else int(loc[0] - each[1] + self.seq_length) for
                              each in rev]
                             for key, rev in self.hmm_locations.items() if not
                             key == hmm]
                reverse_l_flat = list(chain.from_iterable(reverse_l))
                min_rloc.append(min(reverse_l_flat))
            self.locs[hmm].append(min(min_floc))
            self.locs[hmm].append(min(min_rloc))
            #print(self.locs[hmm])
        return self.locs

    def calculate_linkage_scores(self):
        """From dict of with smallest distance locations up- and downstream
        the average is calculated, totaled and a relative value is calculated"""
        try:
            self.locs
        except AttributeError:
            self.find_neighbour_distance()
        linkage_absvals = {hmm: (loc[0] + loc[1]) / 2 for (hmm, loc) in
                           self.locs.items()}
        total_distance = sum([linkVal for hmm, linkVal in
                              linkage_absvals.items()])
        #print(linkage_absvals)
        #print(total_distance)
        linkage_rel_vals = {hmm: [(linkVal / total_distance)]
                            for hmm, linkVal in linkage_absvals.items()}
        #print(self.linkage_rel_vals)
        # Send results to function with lock to write to single file
        # Once all results are there average and make boxplot
        # Write resulting weights to .weights file
        return linkage_rel_vals
