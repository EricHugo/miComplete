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
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqUtils import GC
from operator import itemgetter
from itertools import chain
from termcolor import cprint
import Bio
import argparse
import tempfile
import logging
import logging.handlers
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

def _worker(seqObject, seq_type, argv, q=None, name=None):
    seqObject = ''.join(seqObject)
    base_name = os.path.basename(seqObject).split('.')[0]
    if not name:
        name = base_name
    log_lvl = logging.WARNING
    if argv.debug:
        log_lvl = logging.DEBUG
    elif argv.verbose:
        log_lvl = logging.INFO
    logger = _configure_logger(q, name, log_lvl)
    logger.log(logging.INFO, "Started work on %s" % name)
    if argv.linkage or argv.completeness:
        logger.log(logging.DEBUG, "Creating proteome")
        if re.match("(gb.?.?)|genbank", seq_type):
            logger.log(logging.DEBUG, "gbk-file, will attempt to extract\
                       proteome from gbk")
            proteome = extract_gbk_trans(seqObject)
            # if output file is found to be empty, extract contigs and translate
            if os.stat(proteome).st_size == 0:
                logger.log(logging.DEBUG, "Failed to extract proteome from\
                           gbk, will extract contigs and create proteome\
                           using create_proteome()")
                contigs = get_contigs_gbk(seqObject, name=name)
                proteome = create_proteome(contigs, name)
        elif seq_type == "fna":
            logger.log(logging.DEBUG, "Nucleotide fasta, will translate\
                       using create_proteome()")
            proteome = create_proteome(seqObject, name)
        elif seq_type == "faa":
            logger.log(logging.DEBUG, "Type is amino acid fasta already")
            proteome = seqObject
    else:
        proteome = False
    fastats = parseSeqStats(seqObject, name, seq_type)
    seq_length, all_lengths, GCcontent = fastats.get_length()
    seqstats = (fastats, seq_length, all_lengths, GCcontent)
    if argv.linkage:
        try:
            assert argv.hmms
        except (AssertionError, NameError):
            raise NameError("A set of HMMs must be provided to calculate linkage")
        comp = calcCompleteness(proteome, name, argv.hmms, evalue=argv.evalue,
                                weights=argv.weights, hlist=argv.hlist,
                                linkage=argv.linkage, lenient=argv.lenient,
                                logger=logger)
        hmm_matches, _, total_hmms = comp.get_completeness()
        if argv.hlist:
            comp.print_hmm_lists(directory=argv.hlist)
        try:
            frac_hmm = len(hmm_matches) / len(total_hmms)
        except TypeError:
            frac_hmm = 0
        if frac_hmm < argv.cutoff:
            perc_hmm = frac_hmm * 100
            logger.log(logging.WARNING, "Only %i%% of markers were found in %s,\
                       cannot be used to calculate linkage")
            cprint("Warning:", "red", end=' ', file=sys.stderr)
            print("%i%% of markers were found in %s, cannot be used to calculate\
                   linkage" % (perc_hmm, name), file=sys.stderr)
            return
        linkage = linkageAnalysis(seqObject, name, seq_type, proteome, seqstats,
                                  hmm_matches, argv.debug, q)
        linkage_vals = linkage.calculate_linkage_scores()
        for hmm, match in hmm_matches.items():
            linkage_vals[hmm].append(match)
        q.put(linkage_vals)
        return linkage_vals
    _compile_results(seq_type, name, argv, proteome, seqstats, q, logger)

def _compile_results(seq_type, name, argv, proteome, seqstats, q=None,
                     logger=None):
    """
    Compile results from sequences passed to requested modules and pass
    to Queue for thread safe output. If no Queue given print directly to
    stdout.
    """
    output = []
    if not seq_type == 'faa':
        fastats, seq_length, all_lengths, GC = seqstats
    else:
        fastats, seq_length, all_lengths, GC = "-", "-", "-", "-"
    output.extend((name, seq_length, GC))
    if argv.completeness:
        comp = calcCompleteness(proteome, name, argv.hmms, evalue=argv.evalue,
                                weights=argv.weights, hlist=argv.hlist,
                                linkage=argv.linkage, lenient=argv.lenient,
                                logger=logger)
        filled_hmms, redun_hmms, total_hmms = comp.quantify_completeness()
        try:
            num_hmms = filled_hmms
        except TypeError:
            num_hmms = 0
        output.append(num_hmms)
        try:
            marker_comp = '%0.3f' % (round(num_hmms / total_hmms, 3))
        except ZeroDivisionError:
            marker_comp = 0
        output.append(marker_comp)
        try:
            redundance = '%0.3f' % (round((redun_hmms) / num_hmms, 3))
        except ZeroDivisionError:
            redundance = 0
        output.append(redundance)
        if argv.weights:
            if num_hmms > 0:
                weighted_comp, weighted_redun = comp.attribute_weights(num_hmms)
            else:
                weighted_comp, weighted_redun = 0, 0
            output.append('%0.3f' % weighted_comp)
            output.append('%0.3f' % weighted_redun)
    # only calculate assembly stats if filetype is fna
    if argv.hlist:
        comp.print_hmm_lists(directory=argv.hlist)
    if not re.match("(gb.?.?)|genbank|faa", seq_type):
        N50, L50, N90, L90 = fastats.get_stats(seq_length, all_lengths)
    else:
        N50, L50, N90, L90 = '-', '-', '-', '-'
    output.extend((N50, L50, N90, L90))
    if q:
        q.put(output)
    else:
        if sys.version_info > (3, 0):
            print(*output, sep='\t')
        else:
            print('\t'.join(map(str, output)))

def _configure_logger(q, name, level=logging.WARNING):
    lq = logging.handlers.QueueHandler(q)
    logformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    lq.setFormatter(logformatter)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(lq)
    #logger.log(logging.INFO, "test")
    return logger

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
    all_bias = defaultdict(list)
    while True:
        write_request = q.get()
        cprint(write_request, "green")
        if write_request == 'done':
            break
        if isinstance(write_request, str):
            #logger.handle(write_request)
            continue
        if isinstance(write_request, list):
            if sys.version_info > (3, 0):
                print(*write_request, sep='\t')
            else:
                print('\t'.join(map(str, write_request)))
            continue
        # for testing logging
        continue
        for hmm, match in sorted(write_request.items(), key=lambda e: e[1],
                                 reverse=True):
            {all_bias[hmm].append(1) if float(stats[3]) / float(stats[2]) > 0.1
             else all_bias[hmm].append(0) for
             stats in match[1]}
            weight = (str(hmm) + '\t' + str(match[0]) + '\n')
            weights_tmp.write(weight)
        weights_tmp.write('-\n')
        weights_tmp.flush()
    for hmm, bias in all_bias.items():
        total_fraction_bias = sum(bias) / len(bias)
        if total_fraction_bias > 0.5:
            #logger.log(logging.WARNING, "More than 50% of found marker %s had\
            #           more 10% of score bias. Consider not using this marker"
            #           % hmm)
            cprint("Warning:", "red", file=sys.stderr, end=' ')
            print("More than %s%% of found markers had a higher than %s%% of" %
                  (0.5 * 100, 0.1 * 100), file=sys.stderr, end=' ')
            print("score bias in marker %s. Consider not using this marker." %
                  hmm, file=sys.stderr)
    weights_tmp.close()
    return weights_file

def weights_output(weights_file, logger=None):
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

def create_proteome(fasta, base_name=None, logger=None):
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
    if not base_name:
        base_name = os.path.basename(fasta).split('.')[0]
    prot_filename = base_name + "_prodigal.faa"
    if sys.version_info > (3, 4):
        subprocess.run(['prodigal', '-i', fasta, '-a', prot_filename],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        subprocess.call(['prodigal', '-i', fasta, '-a', prot_filename],
                        stdout=open(os.devnull, 'wb'), stderr=subprocess.STDOUT)
    return prot_filename

def extract_gbk_trans(gbkfile, outfile=None, logger=None):
    """Extract translated sequences from given GeneBank file and out to given
    filename. Returns file name of translation"""
    input_handle = open(gbkfile, mode='r')
    if outfile:
        output_handle = open(outfile, mode='w+')
    else:
        base_name = os.path.basename(gbkfile).split('.')[0]
        outfile = base_name + "_translations.faa"
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
                    #conitnue # not sure why it was here. miscopied from above?
                    output_handle.write(header)
                    output_handle.write(' # ' + start + ' # ' + end + ' # ' +
                                        strand + ' # ID=' + str(contig_n) + '_'
                                        + str(cds_n) + ';')
                    output_handle.write(fasta_trans)
                    continue
                # very occasionally translated seqs have two or more locations
                # regex search to handle such cases
                output_handle.write(header)
                locs = loc_search.search(str(feature.location))
                if locs:
                    locs_list = locs.group(1).split(',')
                    for loc in locs_list:
                        loc_str = ''.join(l for l in ''.join(loc)
                                          if l not in "[]")
                        loc_str = re.sub('\(|\)|:', ' # ', loc_str)
                        loc_str = re.sub('>|<', '', loc_str)
                        output_handle.write(" # " + loc_str)
                else:
                    loc_str = ''.join(l for l in ''.join(str(feature.location))
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

def get_contigs_gbk(gbk, name=None, logger=None):
    """Extracts all sequences from gbk file, returns filename"""
    handle = open(gbk, mode='r')
    if not name:
        name = os.path.basename(gbk).split('.')[0]
    out_handle = open(name, mode='w')
    for seq in SeqIO.parse(handle, "genbank"):
        out_handle.write(">" + seq.id + "\n")
        out_handle.write(str(seq.seq) + "\n")
    out_handle.close()
    return name

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
            action='store_true', help="""Perform completeness check (also requires
            a set of HMMs to have been provided""") 
    parser.add_argument("--lenient", action='store_true', default=False,
            help="""By default miComplete drops hits with too high bias 
            or too low best domain score. This argument disables that behavior, 
            permitting any hit that meets the evalue requirements.""")
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
        input_seqs = [seq.strip().split('\t') for seq in seq_file
                      if not re.match('#|\n', seq)]
    
    if sys.version_info > (3, 0):
        print(*input_seqs, sep='\n', file=sys.stderr)
    else:
        for each_input in input_seqs: print(each_input, file=sys.stderr)

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
    #logfile = logging.FileHandler(args.log, mode='w+')
    logger = _configure_logger(q, "main", "INFO")
    logger.log(logging.INFO, "miComplete has started")
    jobs = []
    for i in input_seqs:
        if len(i) == 2:
            i.append(None)
        job = pool.apply_async(_worker, (i[0], i[1], args),
                               {"q":q, "name":i[2]})
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

if __name__ == "__main__":
    main()
