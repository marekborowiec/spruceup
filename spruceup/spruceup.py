#! /usr/bin/env python3

# coding: utf-8
import configparser
import os
import random
import json
import pdb
from sys import argv, exit
import multiprocessing as mp
from functools import partial
from math import log

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scp
from tqdm import tqdm
import psutil

import aln_parsing, aln_writing

plt.switch_backend('agg')


def read_config(config_file_name):
    try:
        with open(config_file_name) as cf:
            config = configparser.RawConfigParser()
            config.read(config_file_name)
    except IOError as ex:
        exit('Sorry, could not open the file: ' + ex.strerror)
    return config


def distances_wrapper(parsed_alignments, cores, data_type, method='p-distance', fraction=1):
    """Use multiple cores to get p-distances from list of alignment dicts.
    
    Keyword args:
    method (str) -- 'p-distance', 'jc69', or 'missing' (default 'p-distance')"""
    if int(cores) == 1:
        for aln_tuple in tqdm(parsed_alignments, desc='Calculating distances'):
            yield get_distances(aln_tuple, method, fraction, data_type)
    elif int(cores) > 1:
        with mp.Pool(processes=cores) as pool:
            with tqdm(total=len(parsed_alignments)) as pbar:
                for i, output in tqdm(
                    enumerate(
                        pool.imap_unordered(
                            partial(get_distances, method=method, fraction=fraction, data_type=data_type),
                            parsed_alignments,
                        )
                    ),
                    desc='Calculating distances',
                ):
                    pbar.update()
                    yield output


def p_distance(seq1, seq2):
    """Calculate p-distance for two sequences.

    Return the Hamming distance between equal-length sequences
    divided by seq length excluding missing data."""
    if len(seq1) != len(seq2):
        raise ValueError('Sequences are of unequal length. Did you align them?')
    eff_len1 = len(seq1.strip('-').strip('?'))
    eff_len2 = len(seq2.strip('-').strip('?'))
    if eff_len1 != 0 and eff_len2 != 0:
        p_distance = sum(
                el1 != el2
                for el1, el2 in zip(seq1, seq2)
                if el1 is not '-'
                and el2 is not '-'
                and el1 is not '?'
                and el2 is not '?'
                )
    else:
        p_distance = 0
    return p_distance, 5

def get_distances(aln_tuple, method, fraction, data_type):
    """Calculate distances or p-distances for alignment.

    Given tuple (alignment name, alignment distances dict)
    return tuple of (alignment name, list of pairwise distances).
    Keyword args:
    method (str) -- 'p-distance' or 'jc69'
    """

    aln_name, aln_dict = aln_tuple
    seqs_to_compare_to = random.sample(
        aln_dict.items(), int(len(aln_dict.items()) * fraction)
    )

    if method == 'p-distance':
        distances = [
            (sp1, sp2, p_distance(seq1, seq2))
            for sp2, seq2 in seqs_to_compare_to
            for sp1, seq1 in aln_dict.items()
        ]
    elif method == 'jc69':
        distances = [
            (sp1, sp2, jc69_correction(p_distance(seq1, seq2), data_type))
            for sp2, seq2 in seqs_to_compare_to
            for sp1, seq1 in aln_dict.items()
        ]
    return (aln_name, distances)


def jc69_correction(p_distance, data_type):
    """Get Jukes-Cantor corrected distances for nucleotides."""
    if p_distance == 0:
        jc69_corrected = 0
    else:
        if data_type == 'nt':
            jc69_corrected = 3 / 4 * log(1 - 4 / 3 * -p_distance)
        if data_type == 'aa':
            jc69_corrected = 19 / 20 * log(1 - 20 / 19 * -p_distance)
    return jc69_corrected


def dist_taxa_wrapper(dist_tuples):
    """Wrapper for getting aligned lists of taxa and distances.

    Given list of tuples of multiple alignments
    return tuple (alignment name : (taxa rows, distances list))."""
    for aln_name, distances in dist_tuples:
        yield (aln_name, get_dist_and_taxa_lists(distances))


def get_dist_and_taxa_lists(distances):
    """Get aligned lists of taxa and distances.

    Given pairwise distances tuples of 
    (taxon1, taxon2, pairwise distance between the two)
    return tuple of (taxa, distances)."""
    taxa_rows = [sp2 for (sp1, sp2, dist) in distances]
    dists = [dist for (sp1, sp2, dist) in distances]
    return (taxa_rows, dists)


def get_dist_matrix(distances, taxa_no):
    """Create distance matrix for alignment."""
    matrix = [distances[x : x + taxa_no] for x in range(0, len(distances), taxa_no)]
    return matrix


def get_taxon_map(taxa_rows):
    """Create taxon map for distance matrix.

    Taxon map is a list of unique taxon names in order.
    """
    seen = set()
    taxa_map = [sp for sp in taxa_rows if not (sp in seen or seen.add(sp))]
    return taxa_map


def mean_distances_wrapper(aln_tpls):
    """Wrapper to get taxon mean distances across multiple alignments.

    Given list of tuples (alignment name, (taxa rows, distances list))
    return list of tuples (alignment name, {taxon : mean distance within alignment}).
    """
    for aln_tpl_lists in aln_tpls:
        aln_name, dist_tuple = aln_tpl_lists
        taxa_rows, dists = dist_tuple
        taxon_map = get_taxon_map(taxa_rows)
        dist_matrix = get_dist_matrix(dists, len(taxon_map))
        means = get_mean_distances(dist_matrix, taxon_map)
        yield (aln_name, means)


def get_mean_distances(dist_matrix, taxon_map):
    """Get taxon mean distances from taxon map and distance matrix.

    Given matrix of all pairwise distances and taxon map
    return dict of {taxon : mean distances}.
    """
    mean_distances = {
        sp: mean_dist
        for sp, mean_dist in zip(
            taxon_map, [get_list_mean(dist) for dist in dist_matrix]
        )
    }
    return mean_distances


def get_list_mean(lst):
    """Return mean for all items in a list."""
    list_mean = sum(lst) / float(len(lst))
    return round(list_mean, 5)


def dists_per_taxon(means_tuple_list):
    """Get mean pairwise distances for taxon in alignment.

    Given tuple list [(alignment name, {taxon : mean distance within alignment})]
    return dict of {taxon : (alignment name, mean distance within alignment)}.
    """
    taxa_dists = {}
    for aln_name, dist_dict in sorted(means_tuple_list):
        for sp, mean_dist in dist_dict.items():
            if sp not in taxa_dists.keys():
                taxa_dists[sp] = [(aln_name, mean_dist)]
            else:
                taxa_dists[sp].append((aln_name, mean_dist))
    return taxa_dists


def means_per_taxon(taxa_dists):
    """Get mean distances for taxon.

    Given dict of {taxon : (alignment, mean distance within alignment)}
    return dict of {taxon : mean distance across alignments}.
    """
    taxa_means = {
        sp: get_list_mean([dists for aln_name, dists in aln_dists])
        for sp, aln_dists in taxa_dists.items()
    }
    return taxa_means


def plot_taxon_dists(dists, taxon, criterion, cutoff, cutoff_line, fit_line=0):
    fname = '{}-{}{}.png'.format(taxon, cutoff, criterion)
    plt.figure(num=None, figsize=(12, 6), dpi=150, facecolor='w', edgecolor='k')
    plt.xlim(0, np.nanmax(dists))
    plt.hist(dists[~np.isnan(dists)], bins=100, density=True)
    if fit_line is not 0:
        x = np.linspace(0, np.nanmax(dists), 500)
        plt.plot(x, fit_line)
    plt.title(taxon)
    plt.axvline(cutoff_line, color='k', linestyle='dashed', linewidth=1)
    plt.savefig(fname)
    plt.close()


def get_taxon_dists(taxa_dists, taxon):
    """Get dict of distances for each alignment for a taxon."""
    return {(taxon, aln_name): dist for aln_name, dist in taxa_dists[taxon]}


def get_outliers_wrapper(
    all_taxa_dists, window_size, criterion, cutoff, manual_cutoffs
):
    """Wrapper around outlier identification function.

    Given dict of {taxon : alignment_name, mean_distance}
    return dict of tuples {taxon : (taxon_mean_distance, outlier_sequence_ranges) }
    """
    taxa = sorted(all_taxa_dists.keys())
    if criterion == 'lognorm':
        outliers_dict = {
            taxon: get_lognorm_outliers(
                get_taxon_dists(all_taxa_dists, taxon),
                taxon,
                window_size,
                criterion,
                cutoff,
                manual_cutoffs,
            )
            for taxon in taxa
        }
    if criterion == 'mean':
        outliers_dict = {
            taxon: get_mean_outliers(
                get_taxon_dists(all_taxa_dists, taxon),
                taxon,
                window_size,
                criterion,
                cutoff,
                manual_cutoffs,
            )
            for taxon in taxa
        }
    return outliers_dict


def get_window_tuple(tpl, window_size):
    taxon, aln = tpl
    aln_start = aln 
    aln_end = aln_start + window_size
    aln_tpl = (aln_start, aln_end)
    return aln_tpl


def get_lognorm_outliers(
    taxon_dists, taxon, window_size, criterion, cutoff, manual_cutoffs
):
    """Identify outlier windows in a taxon.

    Given dict of _(taxon, aln_name) : dist}
    return tuple of lognormal fit cutoff for taxon 
    and list of ranges in sequence that are outliers.
    """
    from scipy.optimize import curve_fit
    dists = np.asarray(list(taxon_dists.values()))
    dists[dists == 0] = np.nan
    dists = dists[~np.isnan(dists)]
    shape, loc, scale = scp.lognorm.fit(dists, floc=0)
    logn_cutoff = scp.lognorm.ppf(cutoff, shape, loc, scale)
    x = np.linspace(0, np.nanmax(dists), 500)
    logn_fit_line = scp.lognorm.pdf(x, shape, loc, scale)
    outliers = []
    if manual_cutoffs:
        manual_dict = {}
        for group in manual_cutoffs:
            manual_taxon_name, manual_cutoff_value = group
            manual_cutoff = float(manual_cutoff_value)
            manual_dict[manual_taxon_name] = manual_cutoff
        if taxon in manual_dict.keys():
            plot_taxon_dists(
                dists,
                taxon,
                criterion,
                cutoff,
                manual_dict[taxon],
                fit_line=logn_fit_line,
            )
            for tpl, dist in sorted(taxon_dists.items()):
                if dist >= manual_dict[taxon]:
                    outliers.append(get_window_tuple(tpl, window_size))
        else:
            plot_taxon_dists(
                dists, taxon, criterion, cutoff, logn_cutoff, fit_line=logn_fit_line
            )
            for tpl, dist in sorted(taxon_dists.items()):
                if dist >= logn_cutoff:
                    outliers.append(get_window_tuple(tpl, window_size))
    else:
        plot_taxon_dists(
            dists, taxon, criterion, cutoff, logn_cutoff, fit_line=logn_fit_line
        )
        for tpl, dist in sorted(taxon_dists.items()):
            if dist >= logn_cutoff:
                outliers.append(get_window_tuple(tpl, window_size))
    if outliers:
        merged_outliers = merge(outliers)
    else:
        merged_outliers = []
    outlier_sequence_ranges = list(merged_outliers)
    return (logn_cutoff, outlier_sequence_ranges)


def get_mean_outliers(
    taxon_dists, taxon, window_size, criterion, cutoff, manual_cutoffs
):
    """Identify outlier windows in a taxon.

    Given dict of _(taxon, aln_name) : dist}
    return tuple of mean cutoff distance for taxon 
    and list of ranges in sequence that are outliers.
    """
    dists = np.asarray(list(taxon_dists.values()))
    dists[dists == 0] = np.nan
    dists = dists[~np.isnan(dists)]    
    mean = np.mean(list(taxon_dists.values()))
    mean_cutoff = round((mean * cutoff), 5)
    outliers = []
    if manual_cutoffs:
        manual_dict = {}
        for group in manual_cutoffs:
            manual_taxon_name, manual_cutoff_value = group
            manual_cutoff = float(manual_cutoff_value)
            manual_dict[manual_taxon_name] = manual_cutoff
        if taxon in manual_dict.keys():
            plot_taxon_dists(
                dists, taxon, criterion, cutoff, manual_dict[taxon]
            )
            for tpl, dist in sorted(taxon_dists.items()):
                if dist >= manual_dict[taxon]:
                    outliers.append(get_window_tuple(tpl, window_size))
        else:
            plot_taxon_dists(
            dists, taxon, criterion, cutoff, mean_cutoff
            )
            for tpl, dist in sorted(taxon_dists.items()):
                if dist >= mean_cutoff:
                    outliers.append(get_window_tuple(tpl, window_size))
    else:
        plot_taxon_dists(dists, taxon, criterion, cutoff, mean_cutoff)
        for tpl, dist in sorted(taxon_dists.items()):
            if dist >= mean_cutoff:
                outliers.append(get_window_tuple(tpl, window_size))
    if outliers:
        merged_outliers = merge(outliers)
    else:
        merged_outliers = []
    outlier_sequence_ranges = list(merged_outliers)
    return (mean_cutoff, outlier_sequence_ranges)


def merge(ranges):
    """Merge a list of overlapping ranges."""
    saved = list(ranges[0])
    for st, en in sorted([sorted(t) for t in ranges]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
            yield tuple(saved)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en


def get_windows(parsed_alignment, window_size, stride):
    # extract alignment windows of desired length and stride
    print('Splitting into size-{} windows ...\n'.format(window_size))
    aln_len = len(next(iter(parsed_alignment.values())))  # random seq length
    # initiate list of window dicts
    list_of_windows = []
    add_to_list_of_windows = list_of_windows.append
    for i in range(0, aln_len, stride):
        # loop over all parsed partitions, adding taxa and sliced sequences
        start = i
        stop = i + window_size
        new_dict = {}
        if stop < aln_len:
            for taxon, seq in parsed_alignment.items():
                new_seq = '{}'.format(seq[start:stop])
                new_dict[taxon] = new_seq
        else:
            for taxon, seq in parsed_alignment.items():
                new_seq = '{}'.format(seq[start:aln_len])
                new_dict[taxon] = new_seq
            break
        add_to_list_of_windows((i, new_dict))
    return list_of_windows


def replace_seq(text, start, end, replacement=''):
    length = end - start + 1
    return '{}{}{}'.format(text[:start], replacement * length, text[end + 1 :])


def print_mem():
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1000000
    print('\n')
    print('Used {0:.2f} MB memory\n'.format(mem))


def remove_outliers(parsed_alignment, outliers_dict):
    aln_name, aln_dict = parsed_alignment
    total_sites_removed = 0
    if outliers_dict:
        trimmed_aln_dict = {}
        for taxon, seq in sorted(aln_dict.items()):
            cutoff_value, ranges = outliers_dict[taxon]
            if ranges:
                for index, r in enumerate(ranges):
                    start, end = r
                    total_sites_removed += end - start
                    if index == 0:
                        new_seq = replace_seq(seq, start, end, '?')
                    else:
                        new_seq = replace_seq(new_seq, start, end, '?')
            else:
                new_seq = seq
            trimmed_aln_dict[taxon] = new_seq
    else:
        trimmed_aln_dict = aln_dict
        # print a warning message?
    return (total_sites_removed, trimmed_aln_dict)


def get_alignment_size(alignment_tuple):
    alignment_name, alignment_dict = alignment_tuple
    seq_length = len(next(iter(alignment_dict.values())))
    total_alignment_size = seq_length * len(alignment_dict.values())
    return total_alignment_size


def get_removed_fraction(untrimmed_alignment_size, no_sites_trimmed):
    removed_fraction = no_sites_trimmed / untrimmed_alignment_size
    return removed_fraction


def print_report(outliers, criterion, cutoff, manual_cutoffs):
    report_string = ''
    for taxon, tpl in sorted(outliers.items()):
        cutoff_value, outliers_list = tpl
        if manual_cutoffs:
            for group in manual_cutoffs:
                manual_taxon_name, manual_cutoff_value = group
                if taxon == manual_taxon_name:
                    cutoff_value = manual_cutoff_value
        ranges = ''
        total_seq_removed_from_taxon = 0
        for outlier_range in outliers_list:
            start, end = outlier_range
            seq_removed = int(end) - int(start)
            one_range = '{}-{}\t'.format(start, end)
            ranges += one_range
            total_seq_removed_from_taxon += seq_removed
        report_string += (
            '{}:\n'
            'Cutoff: {}\n'
            'Removed {} positions\n'
            '{}\n\n'.format(taxon, cutoff_value, total_seq_removed_from_taxon, ranges)
        )
    return report_string


def write_report(report_string, report_file_name):
    with open(report_file_name, 'w') as rf:
        rf.write(report_string)

def write_distances_dict(mean_taxon_distances, window_size, overlap):
    dist_fn = 'distances-{}window-{}overlap.json'.format(window_size, overlap)
    with open(dist_fn, 'w') as fp:
        print('\n')
        print('Writing distances to file {}'.format(dist_fn))
        print('\n')
        json.dump(mean_taxon_distances, fp)


def read_distances_dict(distances_json):
    with open(distances_json, 'r') as fp:
        print('\n')
        print('Reading distances from file {}'.format(distances_json))
        print('\n')
        mean_taxon_distances = json.load(fp)
        return mean_taxon_distances

def analyze(
    alignment_file_name, input_file_format, data_type, window_size, overlap, cores, method, fraction
):
    print('Parsing alignment {} ...\n'.format(alignment_file_name))
    aln_tuple = aln_parsing.parse_alignment(alignment_file_name, input_file_format)
    stride = get_stride(window_size, overlap)
    aln_name, aln_dict = aln_tuple
    windows = get_windows(aln_dict, window_size, stride)
    all_distances = distances_wrapper(windows, cores, data_type, method=method, fraction=fraction)
    taxa_distances = dist_taxa_wrapper(all_distances)
    mean_aln_distances = mean_distances_wrapper(taxa_distances)
    mean_taxon_distances = dists_per_taxon(mean_aln_distances)
    write_distances_dict(mean_taxon_distances, window_size, overlap) 
    return (aln_tuple, mean_taxon_distances)


def output_loop(
    untrimmed_alignment,
    distances,
    window_size,
    criterion,
    cutoffs,
    manual_cutoffs,
    report_file_name,
    out_format,
    out_file_name,
    data_type,
):
    alignment_sites = get_alignment_size(untrimmed_alignment)
    for cutoff_string in cutoffs:
        cutoff = float(cutoff_string)
        print('Finding outliers for {} {}s cutoff ...'.format(cutoff, criterion))
        outliers = get_outliers_wrapper(
            distances, window_size, criterion, cutoff, manual_cutoffs
        )
        sites_removed, trimmed_alignment = remove_outliers(
            untrimmed_alignment, outliers
        )
        print('Sites removed: {}'.format(sites_removed))
        percent_removed = get_removed_fraction(alignment_sites, sites_removed) * 100
        print(
            'Removed {:.2f}% of sites at cutoff of {} {}s'.format(
                percent_removed, cutoff_string, criterion
            )
        )
        cutoff_report_fname = '{}_{}s-cutoff-{}'.format(
            cutoff, criterion, report_file_name
        )
        cutoff_trimmed_aln_fname = '{}_{}s-cutoff-{}'.format(
            cutoff, criterion, out_file_name
        )
        write_report(
            print_report(outliers, criterion, cutoff, manual_cutoffs),
            cutoff_report_fname,
        )
        print('Wrote report {} ...'.format(cutoff_report_fname))
        aln_writing.write_alignment_file(
            trimmed_alignment, out_format, cutoff_trimmed_aln_fname, data_type
        )
        print('Wrote trimmed alignment {} ...\n'.format(cutoff_trimmed_aln_fname))

def get_stride(window_size, overlap):
    return window_size - overlap

def main():
    script, config_file_name = argv
    conf = read_config(config_file_name)
    # input
    alignment_name = conf.get('input', 'input_file_name')
    file_format = conf.get('input', 'input_format')
    data_type = conf.get('input', 'data_type')
    distances_json = conf.get('input', 'distances_object_file')
    # analysis
    method = conf.get('analysis', 'distance_method')
    window_size = conf.getint('analysis', 'window_size')
    overlap = conf.getint('analysis', 'overlap')
    # include warning if large stride will result in non-overlapping windows
    fraction = conf.getfloat('analysis', 'fraction')
    cores = conf.getint('analysis', 'cores')
    criterion = conf.get('analysis', 'criterion')
    cutoffs = conf.get('analysis', 'cutoffs').split(',')
    manual_cutoffs = conf.get('analysis', 'manual_cutoffs')
    if not manual_cutoffs: # this should be resolved differently
        manual_cutoffs = False
    else:
        manual_cutoffs = [
            tuple(taxon_cutoff.split(','))
            for taxon_cutoff in conf.get('analysis', 'manual_cutoffs').split(';')
        ]
    # output
    output_file_aln = conf.get('output', 'output_file_aln')
    output_format = conf.get('output', 'output_format')
    report = conf.get('output', 'report')
    if distances_json:
        print('Parsing alignment {} ...\n'.format(alignment_name))
        alignment = aln_parsing.parse_alignment(alignment_name, file_format)
        mean_taxon_distances = read_distances_dict(distances_json)
    else:
        alignment, mean_taxon_distances = analyze(
            alignment_name, file_format, data_type, window_size, overlap, cores, method, fraction
        )

    print_mem()
    output_loop(
        alignment,
        mean_taxon_distances,
        window_size,
        criterion,
        cutoffs,
        manual_cutoffs,
        report,
        output_format,
        output_file_aln,
        data_type,
    )


if __name__ == '__main__':

    main()
### To do:

# 1) add docstrings to all functions
# 2) input validation
# 3) tests
# 4) setup.py
