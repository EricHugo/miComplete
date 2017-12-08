#!/usr/bin/env python

"""
Copyright and License

Copyright 2017 Eric Hugoson and Lionel Guy

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. 

See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with 
this program.  If not, see http://www.gnu.org/licenses/.
"""

from __future__ import print_function, division
from distutils import spawn
from collections import defaultdict, namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from operator import itemgetter
from itertools import chain
from termcolor import cprint
import Bio
import argparse
import tempfile
import logging
import sys
import shutil
import subprocess
import math
import re
import os
import threading
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
try:
    from statistics import median
except ImportError:
    from numpy import median
try:
    import queue
except ImportError:
    import Queue as queue
try:
    from micomplete import parseSeqStats
    from micomplete import linkageAnalysis
    from micomplete import calcCompleteness
except ImportError:
    from parseseqs import parseSeqStats
    from linkageanalysis import linkageAnalysis
    from completeness import calcCompleteness

def _worker(seqObject, seqType, argv, q=None, name=None):
    seqObject = ''.join(seqObject)
    baseName = os.path.basename(seqObject).split('.')[0]
    if not name:
        name = baseName
    if argv.linkage or argv.completeness:
        if re.match("(gb.?.?)|genbank", seqType):
            proteome = extract_gbk_trans(seqObject)
        elif seqType == "fna":
            proteome = create_proteome(seqObject, name)
        elif seqType == "faa":
            proteome = seqObject
    else:
        proteome = False
    fastats = parseSeqStats(seqObject, name, seqType)
    seqLength, allLengths, GCcontent = fastats.get_length()
    seqstats = (fastats, seqLength, allLengths, GCcontent)
    if argv.linkage:
        try:
            assert argv.hmms
        except (AssertionError, NameError):
            raise NameError("A set of HMMs must be provided to calculate linkage")
        comp = calcCompleteness(proteome, name, argv.hmms, argv.evalue, 
                argv.weights, argv.hlist, argv.linkage, argv.debug)
        hmmMatches, dupHmms, totalHmms = comp.get_completeness()
        if argv.hlist:
            comp.print_hmm_lists(directory=argv.hlist)
        try:
            fracHmm = len(hmmMatches) / len(totalHmms)
        except TypeError:
            fracHmm = 0
        if fracHmm < argv.cutoff:
            percHmm = fracHmm * 100
            cprint("Warning:", "red", end=' ', file=sys.stderr)
            print("%i%% markers were found in %s, cannot be used to calculate linkage" 
                    % (percHmm, name), file=sys.stderr)
            return
        linkage = linkageAnalysis(seqObject, name, seqType, 
                proteome, seqstats, hmmMatches, argv.debug, q)
        linkageVals = linkage.calculate_linkage_scores()
        q.put(linkageVals)
        return linkageVals
    else:
        compile_results(seqType, name, argv, proteome, seqstats, q)

def compile_results(seqType, name, argv, proteome, seqstats, q=None):
    """
    Compile results from sequences passed to requested modules and pass
    to Queue for thread safe output. If no Queue given print directly to
    stdout.
    """
    output = []
    # unpack tuple
    if not seqType == 'faa':
        fastats, seqLength, allLengths, GC = seqstats
    else:
        fastats, seqLength, allLengths, GC = "-", "-", "-", "-"
    output.extend((name, seqLength, GC))
    if argv.completeness:
        comp = calcCompleteness(proteome, name, argv.hmms, argv.evalue, 
                argv.weights, argv.hlist, argv.linkage, argv.debug)
        filledHmms, redunHmms, totalHmms = comp.quantify_completeness()
        try:
            numHmms = filledHmms
        except TypeError:
            numHmms = 0
        output.append(numHmms)
        try:
            markerComp = '%0.3f' % (round(numHmms / totalHmms, 3))
        except ZeroDivisionError:
            markerComp = 0
        output.append(markerComp)
        try:
            redundance = '%0.3f' % (round((redunHmms) / numHmms, 3))
        except ZeroDivisionError:
            redundance = 0
        output.append(redundance)
        if argv.weights:
            if numHmms > 0:
                weightedComp, weightedRedun = comp.attribute_weights(numHmms)
            else:
                weightedComp, weightedRedun = 0, 0
            output.append('%0.3f' % weightedComp)
            output.append('%0.3f' % weightedRedun)
    # only calculate assembly stats if filetype is fna
    if argv.hlist:
        comp.print_hmm_lists(directory=argv.hlist)
    if not re.match("(gb.?.?)|genbank|faa", seqType): 
        N50, L50, N90, L90 = fastats.get_stats(seqLength, allLengths)
    else:
        N50, L50, N90, L90 = '-', '-', '-', '-'
    output.extend((N50, L50, N90, L90))
    if q:
        q.put(output)
    else:
        if sys.version_info > (3, 0):
            print(*write_request, sep='\t')
        else:
            print('\t'.join(map(str, write_request)))

def _listener(q):
    """
    Function responsible for outputting information in a thread safe manner.
    Recieves write requests from Queue and writes different targets depending
    on type. Writes results from general operation to stdout, logobjects
    to logfile and any calculated weights of each organism to unified 
    tmp file.
    """
    weights_file = "micomplete_weights.temp"
    weights_tmp = open(weights_file, mode='w+')
    while True:
        write_request = q.get()
        #print(write_request)
        if write_request == 'done':
            break
        if type(write_request) is str:
            #logger.handle(write_request)
            continue
        if type(write_request) is list:
            if sys.version_info > (3, 0):
                print(*write_request, sep='\t')
            else:
                print('\t'.join(map(str, write_request)))
            continue
        for hmm, weight in sorted(write_request.items(), 
                key=lambda e: e[1], reverse=True):
            weight = (str(hmm) + '\t' + str(weight) + '\n')
            weights_tmp.write(weight)
        weights_tmp.write('-\n')
        weights_tmp.flush()
        boxplot = True
    weights_tmp.close()
    return weights_file

def weights_output(weights_file):
    """Creates boxplot all linkage values for each marker present in
    micomplete_weights.temp file"""
    # also output weights
    # weight = median / sum(medians)
    with open(weights_file, mode='r') as weights:
        hmm_weights = defaultdict(list)
        for weight in weights:
            if weight.strip() == '-':
                continue
            weight = weight.split()
            hmm_weights[weight[0]].append(float(weight[1]))
    data = []
    labels = []
    median_weights = {}
    # establish medians, also append axis data
    for hmm, weight in hmm_weights.items():
        median_weights[hmm] = median(weight)
        data.append(weight)
        labels.append(hmm)
    # calculate normalized median weights
    weights_sum = sum(median_weights.values())
    for hmm, median_weight in median_weights.items():
        norm_weight = median_weight / weights_sum
        print(hmm + "\t" + str(norm_weight))
    # create boxplot
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    plt.boxplot(data, vert=False, sym='+')
    ax.set_yticklabels(labels, fontsize=8)
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
            alpha=0.5)
    ax.set_xlabel('Relative linkage distance')
    ax.set_ylabel('Markers')
    plt.show()

def create_proteome(fasta, baseName=None):
    """Create proteome from given .fna file, returns proteome filename"""
    try:
        assert shutil.which('prodigal')
    except AssertionError:
        raise RuntimeError("Unable to locate prodigal in path")
    except AttributeError:
        try:
            assert spawn.find_executable('prodigal')
        except AssertionError:
            raise RuntimeError('Unable to find prodigal in path')
    if not baseName:
        baseName = os.path.basename(fasta).split('.')[0]
    protFileName = baseName + "_prodigal.faa"
    if sys.version_info > (3, 4):
        subprocess.run(['prodigal', '-i', fasta, '-a', protFileName], 
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        subprocess.call(['prodigal', '-i', fasta, '-a', protFileName], 
                stdout=open(os.devnull, 'wb'), stderr=subprocess.STDOUT)
    return protFileName

def extract_gbk_trans(gbkfile, outfile=None):
    """Extract translated sequences from given GeneBank file and out to given
    filename. Returns file name of translation"""
    input_handle = open(gbkfile, mode='r')
    if outfile:
        output_handle = open(outfile, mode='w+')
    else:
        baseName = os.path.basename(gbkfile).split('.')[0]
        outfile = baseName + "_translations.faa"
        output_handle = open(outfile, mode='w+')
    contig_n = 0
    cds_n = 0
    # compile regex to match multi-loc hits
    loc_search = re.compile('\{(.*?)\}')
    for record in SeqIO.parse(input_handle, "genbank"):
        for feature in record.features:
            if feature.type == "source":
                contig_n += 1
            if feature.type == "CDS":
                cds_n += 1
                try:
                    header = ">" + feature.qualifiers['locus_tag'][0]
                except KeyError:
                    continue
                # some CDS do not have translations, retrieve nucleotide sequence 
                # and translate
                try:
                    assert feature.qualifiers['translation'][0]
                except (AssertionError, KeyError):
                    start = str(feature.location.nofuzzy_start)
                    end = str(feature.location.nofuzzy_end)
                    if feature.strand > 0:
                        strand = '+'
                    else:
                        strand = '-'
                    ttable = int(''.join(feature.qualifiers['transl_table']))
                    # check if the locus represents valid CDS else skip
                    try:
                        fasta_trans = ('\n' + str(feature.extract(record.seq).translate(
                            table=ttable, cds=True, to_stop=True)) + '\n')
                    except Bio.Data.CodonTable.TranslationError:
                        continue
                    continue
                    output_handle.write(header)
                    output_handle.write(' # ' + start + ' # ' + end + ' # ' + 
                            strand + ' # ID=' + str(contig_n) + '_' + str(cds_n)
                            + ';')
                    output_handle.write(fasta_trans)
                # very occasionally translated seqs have two or more locations
                # regex search to handle such cases
                output_handle.write(header)
                locs = loc_search.search(str(feature.location))
                if locs:
                    locs_list = locs.group(1).split(',')
                    for loc in locs_list:
                        loc_str = ''.join( l for l in ''.join(loc) 
                                if l not in "[]")
                        loc_str = re.sub('\(|\)|:', ' # ', loc_str)
                        loc_str = re.sub('>|<', '', loc_str)
                        output_handle.write(" # " + loc_str)
                else:
                    loc_str = ''.join( l for l in ''.join(str(feature.location)) 
                            if l not in "[]")
                    loc_str = re.sub('\(|\)|\:', ' # ', loc_str)
                    loc_str = re.sub('>|<', '', loc_str)
                    output_handle.write(" # " + loc_str)
                output_handle.write('ID=' + str(contig_n) + '_' + str(cds_n) 
                        + ';')
                output_handle.write('\n' + feature.qualifiers['translation'][0]
                        + '\n')
    input_handle.close()
    output_handle.flush()
    output_handle.close()
    return outfile

def main():
    parser = argparse.ArgumentParser(
        description="""
            Quality control of metagenome assembled genomes. Able to gather
            relevant statistics of completeness and redundance of genomes given 
            sets of marker genes, as well as weighted versions of these statstics
            (including defining new weights for any given set).
            """,
        epilog="""Report issues and bugs to the issue tracker at 
                https://bitbucket.org/evolegiolab/micomplete or directly to 
                eric@hugoson.org""")

    parser.add_argument("sequence", help="""Sequence(s) along with type (fna, 
            faa, gbk) provided in a tabular format""")
    parser.add_argument("-t", "--total", required=False, default=False,
            action='store_true', help="""Print total (not implemented)""")
    parser.add_argument("-c", "--completeness", required=False, default=False,
            action='store_true', help="""Do completeness check (also requires
            a set of HMMs to have been provided""") 
    parser.add_argument("--hlist", required=False, default=None, type=str,
            nargs='?', help="""Write list of Present, Absent and
            Duplicated markers for each organism to file""")
    parser.add_argument("--hmms", required=False, default=False,
            help="""Specifies a set of HMMs to be used for completeness check 
            or linkage analysis""")
    parser.add_argument("--weights", required=False, default=False,
            help="""Specify a set of weights for the HMMs specified,
            (optional)""")
    parser.add_argument("--linkage", required=False, default=False, 
            action='store_true', help="""Specifies that the provided sequences 
            should be used to calculate the weights of the provided HMMs""")
    parser.add_argument("--evalue", required=False, type=float, default=1e-10,
            help="""Specify e-value cutoff to be used for completeness check. 
            Default = 1e-10""")
    parser.add_argument("--cutoff", required=False, type=float, default=0.9,
            help="""Specify cutoff percentage of markers required to be present 
            in genome for it be included in linkage calculation. 
            Default = 0.9""")
    parser.add_argument("--threads", required=False, default=1, type=int,
            help="""Specify number of threads to be used in parallel""")
    parser.add_argument("--log", required=False, default="miComplete.log",
            type=str, help="""Log name (default=miComplete.log)""")
    parser.add_argument("-v", "--verbose", required=False, default=False,
            action='store_true', help="""Enable verbose logging""")
    parser.add_argument("--debug", required=False, default=False, 
            action='store_true')
        
    args = parser.parse_args()
    logger = logging.getLogger("miComplete")
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    logfile = logging.FileHandler(args.log, mode='w+')
    logformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logfile.setFormatter(logformatter)
    logger.addHandler(logfile)
    logger.info("miComplete has started")

    if args.completeness or args.linkage:
        try:
            assert shutil.which('hmmsearch')
        except AssertionError:
            raise RuntimeError('Unable to find hmmsearch in path')
        except AttributeError:
            try:
                assert spawn.find_executable('hmmsearch')
            except AssertionError:
                raise RuntimeError('Unable to find hmmsearch in path')

    with open(args.sequence) as seq_file:
        inputSeqs = [ seq.strip().split('\t') for seq in seq_file 
                    if not re.match('#|\n', seq) ]
    
    if sys.version_info > (3, 0):
        print(*inputSeqs, sep='\n', file=sys.stderr)
    else:
        for eachInput in inputSeqs: print(eachInput, file=sys.stderr)

    ## print out column headers, unless linkage is requested
    if args.completeness and not args.linkage:
        if not args.hmms:
            raise TypeError("""Completeness check requires a set of HMMs to 
            look for.""")
        if args.weights:
            print("Name\tLength\tGC-content\tPresent Markers\tCompleteness"
                    "\tRedundance\tCompletenessW\tRedudanceW\tN50\tL50\tN90\tL90")
        else:
            print("Name\tLength\tGC-content\tPresent Markers\tCompleteness"
                    "\tRedundance\tN50\tL50\tN90\tL90")
    elif not args.linkage:
        print("Name\tLength\tGC-content\tN50\tL50\tN90\tL90")

    if mp.cpu_count() < args.threads:
        raise RuntimeError('Specified number of threads are larger than the '
                'number detected in the system: ' + mp.cpu_count())
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads + 1)
    writer = pool.apply_async(_listener, (q,))
    
    jobs = []
    for i in inputSeqs: 
        if len(i) == 2:
            i.append(None)
        job = pool.apply_async(_worker, (i[0], i[1], args, q, i[2]))
        jobs.append(job)

    # get() all processes to catch errors
    for job in jobs:
        job.get()
    q.put("done")
    weights_file = writer.get()
    pool.close()
    pool.join()
    if args.linkage:
        weights_output(weights_file)
    logger.info("miComplete has finished")    

if __name__=="__main__":
    main()
