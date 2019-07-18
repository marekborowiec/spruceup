#! /usr/bin/env python3

# coding: utf-8
import configparser
import logging
import os
import random
import json
import pdb
from sys import argv, exit
import multiprocessing as mp
import time
from functools import partial
from math import log

import treeswift
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as scp
from tqdm import tqdm
import psutil

import aln_parsing, aln_writing

plt.switch_backend('agg')


def get_tree_dist_dict(tree_fn):
    """Make a dict of all-by-all distances from input guide tree."""
    t = treeswift.read_tree_newick(tree_fn)
    taxa_nodes = t.label_to_node()
    tree_dist_dict = {}
    for sp1, node1 in taxa_nodes.items():
        for sp2, node2 in taxa_nodes.items():
            if sp1 != sp2:
                tree_dist = get_tree_dist_between_two_leaves(t, node1, node2)
            else:
                tree_dist = 0
            tree_dist_dict[(sp1, sp2)] = tree_dist
    to_remove = []
    for tpl, distance in tree_dist_dict.items():
        sp1, sp2 = tpl
        if sp1 == '' or sp2 == '':
            to_remove.append(tpl)
    for tpl in to_remove:
        tree_dist_dict.pop(tpl)
    return tree_dist_dict


def get_tree_dist_between_two_leaves(tree, nodeA, nodeB):
    """Given ete phylogenetic tree, get distance between leaves."""
    return tree.distance_between(nodeA, nodeB)


def lookup_tree_dist(tree_dist_dict, sp1, sp2):
    """Fast lookup of tree distance from dict of all guide tree distances."""
    return tree_dist_dict[(sp1, sp2)]


def replace_missing_in_dict(parsed_aln_dict, data_type):
    """Convert ambiguous and missing data to '?' before calculating distances."""
    nt_missing_ambiguous_chars = [
        'K',
        'M',
        'R',
        'Y',
        'S',
        'W',
        'B',
        'V',
        'H',
        'D',
        'X',
        'N',
        'O',
    ]
    aa_missing_ambiguous_chars = ['B', 'J', 'Z', 'X', '.', '*']
    if data_type == 'aa':
        new_dict = {
            taxon: replace_missing_ambiguous(seq, aa_missing_ambiguous_chars)
            for taxon, seq in parsed_aln_dict.items()
        }
    elif data_type == 'nt':
        new_dict = {
            taxon: replace_missing_ambiguous(seq, nt_missing_ambiguous_chars)
            for taxon, seq in parsed_aln_dict.items()
        }
    return new_dict


def replace_missing_ambiguous(seq, missing_ambiguous_list):
    """Given sequence and list of missing or ambiguous characters,
    replace them in sequence with '?'.
    """
    for char in missing_ambiguous_list:
        seq = seq.replace(char, '?')
    return seq


def read_config(config_file_name):
    """Read in configuration file if it exists."""
    try:
        with open(config_file_name) as cf:
            config = configparser.RawConfigParser()
            config.read(config_file_name)
    except IOError as ex:
        exit(
            'Sorry, could not open the configuration file "{}": {}'.format(
                config_file_name, ex.strerror
            )
        )
    return config


def p_distance(seq1, seq2):
    """Calculate Hamming/p-distance for two sequences.

    Return the Hamming distance between equal-length sequences
    or return string NaN if comparing to empty sequence.
    """
    if len(seq1) != len(seq2):
        raise ValueError('Sequences are of unequal length. Did you align them?')
    eff_len1 = len(seq1.strip('-').strip('?'))
    eff_len2 = len(seq2.strip('-').strip('?'))
    if eff_len1 != 0 and eff_len2 != 0:
        p_distance = sum(
            el1 != el2
            for el1, el2 in zip(seq1, seq2)
            if el1 != '-' and el2 != '-' and el1 != '?' and el2 != '?'
        )
    else:
        p_distance = 'NaN'
    return (eff_len1, p_distance)


def get_scaled_distance(distance_tpl):
    """Scale distances from 0 to 1.

    Given tuple of (seq length, Hamming distance)
    return proportion of different sites
    unless one of the sequences was missing.
    """
    eff_seq_len, distance = distance_tpl
    if distance != 'NaN':
        if eff_seq_len > 0 and distance > 0:
            scaled_distance = distance / eff_seq_len
        else:
            scaled_distance = 0
    else:
        scaled_distance = 'NaN'
    return scaled_distance


def jc_correction(distance_tpl, data_type):
    """Get Jukes-Cantor corrected distances.

    Given distance tuple and depending on data type
    compute JC-corrected distance for DNA or proteins and
    return tuple of (seq length, corrected distance).
    """
    eff_seq_len, p_distance = distance_tpl
    if p_distance != 'NaN':
        if data_type == 'nt':
            jc_corrected = 3 / 4 * log(1 - 4 / 3 * -p_distance)
        elif data_type == 'aa':
            jc_corrected = 19 / 20 * log(1 - 20 / 19 * -p_distance)
    else:
        jc_corrected = 'NaN'
    return (eff_seq_len, jc_corrected)


def get_distances_scaled_by_tree(
    method, data_type, tree_dists, sp1, sp2, seq1, seq2
):
    """Given a tree, scale computed distances by tree distances.

    This is done by dividing each spruceup distance 
    by distance between same OTUs on the guide tree.
    """
    if method == 'uncorrected':
        scaled_distance = get_scaled_distance(p_distance(seq1, seq2))
    elif method == 'jc':
        scaled_distance = get_scaled_distance(
            jc_correction(p_distance(seq1, seq2), data_type)
        )
    tree_distance = lookup_tree_dist(tree_dists, sp1, sp2)
    try:
        if scaled_distance == 'NaN':
            tree_scaled = 'NaN'
        else:
            tree_scaled = scaled_distance / tree_distance
    except ZeroDivisionError:
        tree_scaled = scaled_distance
    return tree_scaled


def get_distances_scaled(method, data_type, seq1, seq2):
    """Convert raw to scaled distances.

    Given raw distance between two sequences, 
    scale them from 0 to 1 and in case of Jukes-Cantor distances also
    depending on whether the sequences are nucleotides or amino acids.
    """
    if method == 'uncorrected':
        scaled_distance = get_scaled_distance(p_distance(seq1, seq2))
    elif method == 'jc':
        scaled_distance = get_scaled_distance(
            jc_correction(p_distance(seq1, seq2), data_type)
        )
    return scaled_distance


def get_distances(aln_tuple, tree_dists, method, fraction, data_type):
    """Calculate uncorrected or JC-corrected distances for alignment.

    Given tuple (alignment name, alignment distances dict),
    return tuple of (alignment name, list of pairwise distances).
    """
    aln_name, aln_dict = aln_tuple
    seqs_to_compare_to = random.sample(
        aln_dict.items(), int(len(aln_dict.items()) * fraction)
    )
    if tree_dists == None:
        distances = [
            (sp1, sp2, get_distances_scaled(method, data_type, seq1, seq2))
            for sp2, seq2 in seqs_to_compare_to
            for sp1, seq1 in aln_dict.items()
        ]
    elif tree_dists != None:
        distances = [
            (
                sp1,
                sp2,
                get_distances_scaled_by_tree(
                    method, data_type, tree_dists, sp1, sp2, seq1, seq2
                ),
            )
            for sp2, seq2 in seqs_to_compare_to
            for sp1, seq1 in aln_dict.items()
        ]
    return (aln_name, distances)


def distances_wrapper(
    parsed_alignments, tree_dists, cores, data_type, method, fraction
):
    """Use multiple cores to get distances from list of alignment dicts."""
    if int(cores) == 1:
        for aln_tuple in tqdm(parsed_alignments, desc='Calculating distances'):
            yield get_distances(
                aln_tuple, tree_dists, method, fraction, data_type
            )
    elif int(cores) > 1:
        with mp.Pool(processes=cores) as pool:
            with tqdm(total=len(parsed_alignments)) as pbar:
                for i, output in tqdm(
                    enumerate(
                        pool.imap_unordered(
                            partial(
                                get_distances,
                                tree_dists=tree_dists,
                                method=method,
                                fraction=fraction,
                                data_type=data_type,
                            ),
                            parsed_alignments,
                        )
                    ),
                    desc='Calculating distances',
                ):
                    pbar.update()
                    yield output


def get_dist_and_taxa_lists(distances):
    """Get aligned lists of taxa and distances.

    Given pairwise distances tuples 
    (taxon1, taxon2, pairwise distance between the two)
    return tuple (taxa list, distances list).
    """
    taxa_rows = [sp2 for (sp1, sp2, dist) in distances]
    dists = [dist for (sp1, sp2, dist) in distances]
    return (taxa_rows, dists)


def dist_taxa_wrapper(dist_tuples):
    """Wrapper for getting aligned lists of taxa and distances.

    Given list of tuples of multiple alignments
    return tuple (alignment name : (taxa rows, distances list)).
    """
    for aln_name, distances in dist_tuples:
        yield (aln_name, get_dist_and_taxa_lists(distances))


def get_dist_matrix(distances, taxa_no):
    """Create distance matrix for alignment.

    Split all distances list [sp1 to sp1 dist, sp1 to sp2 dist, sp1 to sp3 dist etc.]
    to return list of lists [[sp1 dists to all other], [sp2 dists to all other], etc.]
    """
    matrix = [
        distances[x : x + taxa_no] for x in range(0, len(distances), taxa_no)
    ]
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
    clean_list = [i for i in lst if i != 'NaN']
    try:
        list_mean = abs(
            sum([i for i in clean_list]) / float(len(clean_list) - 1)
        )  # - 1 ensures that distance to self does not count
    except ZeroDivisionError:  # when there is only one non-empty sequence in window
        list_mean = 0
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


def get_np_dists(dist_list):
    """Get numpy distances array in which zeros are NaN."""
    dists = np.asarray(dist_list)
    dists[dists == 0] = np.nan
    return dists[~np.isnan(dists)]


def get_shape_loc_scale(dists):
    """Given distances fit lognormal distribution."""
    return scp.lognorm.fit(dists, floc=0)


def get_lognorm_fit_line(dists, shape, loc, scale):
    """Get parameters for lognorm distribution plotting."""
    x = np.linspace(0, np.nanmax(dists), 500)
    return scp.lognorm.pdf(x, shape, loc, scale)


def get_lognorm_cutoff(cutoff, shape, loc, scale):
    """Get cutoff values for plotting."""
    logn_cutoff = scp.lognorm.ppf(cutoff, shape, loc, scale)
    return logn_cutoff


def get_mean_cutoff(dist_list, cutoff):
    """Get mean from a list of distances."""
    mean = np.mean(dist_list)
    return round((mean * cutoff), 5)


def plotting_wrapper(
    all_taxa_dists, window_size, method, criterion, cutoffs, manual_cutoffs
):
    """This is a wrapper for plot_taxon_dists() function to work across all OTUs and windows."""
    taxa = sorted(all_taxa_dists.keys())
    for taxon in taxa:
        if criterion == 'lognorm':
            dist_list = [window[1] for window in all_taxa_dists[taxon]]
            dists = get_np_dists(dist_list)
            shape, loc, scale = get_shape_loc_scale(dists)
            logn_fit_line = get_lognorm_fit_line(dists, shape, loc, scale)
            plot_taxon_dists(
                all_taxa_dists,
                taxon,
                method,
                criterion,
                cutoffs,
                fit_line=logn_fit_line,
            )
        if criterion == 'mean':
            plot_taxon_dists(all_taxa_dists, taxon, method, criterion, cutoffs)


def plot_taxon_dists(
    all_taxa_dists, taxon, method, criterion, cutoffs, fit_line=0
):
    """Get a histogram plot of distance distribution across windows."""
    fname = '{}-{}-{}.png'.format(taxon, method, criterion)
    dist_list = [window[1] for window in all_taxa_dists[taxon]]
    dists = get_np_dists(dist_list)
    plt.figure(num=None, figsize=(12, 6), dpi=150, facecolor='w', edgecolor='k')
    plt.xlim(0, np.nanmax(dists))
    plt.hist(dists[~np.isnan(dists)], bins=100, density=True)
    if fit_line is not 0:
        x = np.linspace(0, np.nanmax(dists), 500)
        plt.plot(x, fit_line)
    plt.title(taxon)
    colors = iter(plt.cm.rainbow(np.linspace(0, 1, len(cutoffs))))
    for cutoff in cutoffs:
        if criterion == 'lognorm':
            shape, loc, scale = get_shape_loc_scale(dists)
            cutoff_line = get_lognorm_cutoff(float(cutoff), shape, loc, scale)
        if criterion == 'mean':
            cutoff_line = get_mean_cutoff(dists, float(cutoff))
        color = next(colors)
        plt.axvline(
            cutoff_line,
            color=color,
            label=str(cutoff),
            linestyle='dashed',
            linewidth=1,
        )
    plt.legend(loc='upper right')
    plt.savefig(fname)
    plt.close()


def get_taxon_dists(taxa_dists, taxon):
    """Get dict of distances for each alignment for a taxon."""
    return {(taxon, aln_name): dist for aln_name, dist in taxa_dists[taxon]}


def get_outliers_wrapper(
    all_taxa_dists, window_size, method, criterion, cutoff, manual_cutoffs
):
    """Wrapper around outlier identification function.

    Given dict of {taxon : alignment_name, mean_distance}
    return dict of tuples {taxon : (taxon_mean_distance, outlier_sequence_ranges)}
    """
    taxa = sorted(all_taxa_dists.keys())
    if criterion == 'lognorm':
        outliers_dict = {
            taxon: get_lognorm_outliers(
                all_taxa_dists,
                taxon,
                window_size,
                method,
                criterion,
                cutoff,
                manual_cutoffs,
            )
            for taxon in taxa
        }
    if criterion == 'mean':
        outliers_dict = {
            taxon: get_mean_outliers(
                all_taxa_dists,
                taxon,
                window_size,
                method,
                criterion,
                cutoff,
                manual_cutoffs,
            )
            for taxon in taxa
        }
    return outliers_dict


def get_window_tuple(tpl, window_size):
    """Given a tuple (window_name, window_size) parse start and end to that window sequence."""
    aln, distance = tpl
    aln_start = aln
    aln_end = aln_start + window_size
    aln_tpl = (aln_start, aln_end)
    return aln_tpl


def get_outliers_list(window_dist_list, cutoff):
    """List comprehension to get all windows above certain threshold."""
    return [window for window in window_dist_list if window[1] >= cutoff]


def get_lognorm_outliers(
    all_dists, taxon, window_size, method, criterion, cutoff, manual_cutoffs
):
    """Identify outlier windows in a taxon.

    Given dict of {(taxon, aln_name) : dist}
    return tuple of lognormal fit cutoff for taxon 
    and list of ranges in sequence that are outliers.
    """
    dist_list = [window[1] for window in all_dists[taxon]]
    dists = get_np_dists(dist_list)
    shape, loc, scale = get_shape_loc_scale(dists)
    logn_cutoff = get_lognorm_cutoff(cutoff, shape, loc, scale)
    if manual_cutoffs:
        manual_dict = {}
        for group in manual_cutoffs:
            manual_taxon_name, manual_cutoff_value = group
            manual_cutoff = float(manual_cutoff_value)
            manual_dict[manual_taxon_name] = manual_cutoff
        if taxon in manual_dict.keys():
            outliers_list = sorted(
                get_outliers_list(all_dists[taxon], manual_dict[taxon])
            )
            outliers = [
                get_window_tuple(window, window_size)
                for window in outliers_list
            ]
        else:
            outliers_list = sorted(
                get_outliers_list(all_dists[taxon], logn_cutoff)
            )
            outliers = [
                get_window_tuple(window, window_size)
                for window in outliers_list
            ]
    else:
        outliers_list = sorted(get_outliers_list(all_dists[taxon], logn_cutoff))
        outliers = [
            get_window_tuple(window, window_size) for window in outliers_list
        ]
    if outliers:
        merged_outliers = merge(outliers)
    else:
        merged_outliers = []
    outlier_sequence_ranges = list(merged_outliers)
    return (logn_cutoff, outlier_sequence_ranges)


def get_mean_outliers(
    all_dists, taxon, window_size, method, criterion, cutoff, manual_cutoffs
):
    """Identify outlier windows in a taxon.

    Given dict of {(taxon, aln_name) : dist}
    return tuple of mean cutoff distance for taxon 
    and list of ranges in sequence that are outliers.
    """
    dist_list = [window[1] for window in all_dists[taxon]]
    dists = get_np_dists(dist_list)
    mean_cutoff = get_mean_cutoff(dists, cutoff)
    if manual_cutoffs:
        manual_dict = {}
        for group in manual_cutoffs:
            manual_taxon_name, manual_cutoff_value = group
            manual_cutoff = float(manual_cutoff_value)
            manual_dict[manual_taxon_name] = manual_cutoff
        if taxon in manual_dict.keys():
            outliers_list = sorted(
                get_outliers_list(all_dists[taxon], manual_dict[taxon])
            )
            outliers = [
                get_window_tuple(window, window_size)
                for window in outliers_list
            ]
        else:
            outliers_list = sorted(
                get_outliers_list(all_dists[taxon], mean_cutoff)
            )
            outliers = [
                get_window_tuple(window, window_size)
                for window in outliers_list
            ]
    else:
        outliers_list = sorted(get_outliers_list(all_dists[taxon], mean_cutoff))
        outliers = [
            get_window_tuple(window, window_size) for window in outliers_list
        ]
    if outliers:
        merged_outliers = merge(outliers)
    else:
        merged_outliers = []
    outlier_sequence_ranges = list(merged_outliers)
    return (mean_cutoff, outlier_sequence_ranges)


def merge(ranges):
    """Merge a list of overlapping ranges.

    e.g. given [[0,20], [5,25], [30,50]] 
    return [[0,25], [30,50]].
    """
    merged = []
    for higher in ranges:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = (lower[0], upper_bound)
            else:
                merged.append(higher)
    return merged


def get_windows(parsed_alignment, window_size, overlap):
    """Split alignment into sequence windows.

    Given dict {taxon: sequence} return list of tuples (window_name, {taxon: window_seq}).
    """
    # extract alignment windows of desired length and stride
    logging.info(
        'Splitting into size-{} windows with {} overlap ...\n'.format(
            window_size, overlap
        )
    )
    aln_len = len(next(iter(parsed_alignment.values())))  # random seq length
    if window_size > aln_len:
        exit(
            'Invalid window size: "{}" is greater than your alignment length ({}).'.format(
                window_size, aln_len
            )
        )
    stride = get_stride(window_size, overlap)
    aln_len_window = aln_len + window_size  # for iteration
    # initiate list of window dicts
    list_of_windows = []
    add_to_list_of_windows = list_of_windows.append
    for i in range(0, aln_len_window, stride):
        # loop over all parsed partitions, adding taxa and sliced sequences
        start = i
        stop = i + window_size
        new_dict = {}
        if stop <= aln_len:
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
    """Replace slice of string given coordinates and replacement character."""
    length = end - start
    return '{}{}{}'.format(text[:start], replacement * length, text[end:])


def print_mem():
    """Print memory usage using psutil module."""
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1_000_000
    logging.info('Used {0:.2f} MB memory\n'.format(mem))


def remove_outliers(parsed_alignment, outliers_dict):
    """Remove outlier sequences and return trimmed alignment.

    Given parsed alignmend dict {taxon : sequence} 
    and outliers dict {taxon : (taxon_mean_distance, outlier_sequence_ranges)}
    trim out sequence identified as outlier 
    and return tuple with count of removed sites and dict {taxon : trimmed_sequence}.
    """
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
    return (total_sites_removed, trimmed_aln_dict)


def get_alignment_size(alignment_tuple):
    """Get alignment length from alignment tuple.

    Given tuple (aln_name, aln_dict) get len of random sequence in {taxon: sequence} alignment dict."""
    alignment_name, alignment_dict = alignment_tuple
    seq_length = len(next(iter(alignment_dict.values())))
    total_alignment_size = seq_length * len(alignment_dict.values())
    return total_alignment_size


def get_removed_fraction(untrimmed_alignment_size, no_sites_trimmed):
    """Calculate percentage of sites removed given alignment size and number of sites trimmed."""
    removed_fraction = no_sites_trimmed / untrimmed_alignment_size
    return removed_fraction


def print_report(outliers, criterion, cutoff, manual_cutoffs):
    """Report per-taxon cutoff value and no of removed alignment positions.

    Construct report string given outliers dict 
    {taxon : (taxon_mean_distance, outlier_sequence_ranges)}."""
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
            '{}\n\n'.format(
                taxon, cutoff_value, total_seq_removed_from_taxon, ranges
            )
        )
    return report_string


def write_report(report_string, report_file_name):
    """Write report string to file."""
    with open(report_file_name, 'w') as rf:
        rf.write(report_string)


def write_distances_dict(
    mean_taxon_distances, distances_method, window_size, overlap
):
    """Write json file with per-taxon windows and their distances.

    The format is: {"taxon": [[window0, mean_distance_in_window], [window1, dist] ...]}.
    """
    dist_fn = '{}-distances-{}window-{}overlap.json'.format(
        distances_method, window_size, overlap
    )
    with open(dist_fn, 'w') as f:
        logging.info('Writing distances to file {} ...\n'.format(dist_fn))
        json.dump(mean_taxon_distances, f)


def read_distances_dict(distances_json):
    """Parse json file with distances."""
    with open(distances_json, 'r') as f:
        logging.info(
            'Reading distances from file {} ...\n'.format(distances_json)
        )
        mean_taxon_distances = json.load(f)
    return mean_taxon_distances


def analyze(
    alignment_file_name,
    tree_file_name,
    input_file_format,
    data_type,
    window_size,
    overlap,
    cores,
    method,
    fraction,
):
    """Load, parse, and analyze the data."""
    logging.info('Parsing alignment {} ...\n'.format(alignment_file_name))
    aln_tuple = aln_parsing.parse_alignment(
        alignment_file_name, input_file_format
    )
    aln_name, aln_dict = aln_tuple
    no_missing_ambiguous_dict = replace_missing_in_dict(aln_dict, data_type)
    if tree_file_name == None:
        logging.info('Did not find guide tree, continuing ...\n')
        tree_dists = None
    else:
        logging.info('Reading in guide tree {} ...\n'.format(tree_file_name))
        tree_dists = get_tree_dist_dict(tree_file_name)
    windows = get_windows(no_missing_ambiguous_dict, window_size, overlap)
    all_distances = distances_wrapper(
        windows, tree_dists, cores, data_type, method=method, fraction=fraction
    )
    taxa_distances = dist_taxa_wrapper(all_distances)
    mean_aln_distances = mean_distances_wrapper(taxa_distances)
    mean_taxon_distances = dists_per_taxon(mean_aln_distances)
    write_distances_dict(mean_taxon_distances, method, window_size, overlap)
    return (aln_tuple, mean_taxon_distances)


def output_loop(
    untrimmed_alignment,
    distances,
    window_size,
    method,
    criterion,
    cutoffs,
    manual_cutoffs,
    report_file_name,
    out_format,
    out_file_name,
    data_type,
):
    """Write and plot output files from analysis."""
    alignment_sites = get_alignment_size(untrimmed_alignment)
    cutoff_floats = [float(cutoff_string) for cutoff_string in cutoffs]
    for cutoff in cutoff_floats:
        logging.info(
            'Finding outliers for {} {}s cutoff ...'.format(cutoff, criterion)
        )
        outliers = get_outliers_wrapper(
            distances, window_size, method, criterion, cutoff, manual_cutoffs
        )
        sites_removed, trimmed_alignment = remove_outliers(
            untrimmed_alignment, outliers
        )
        logging.info('Sites removed: {}'.format(sites_removed))
        percent_removed = (
            get_removed_fraction(alignment_sites, sites_removed) * 100
        )
        logging.info(
            'Removed {:.2f}% of sites at cutoff of {} {}s'.format(
                percent_removed, cutoff, criterion
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
        logging.info('Wrote report {}'.format(cutoff_report_fname))
        aln_writing.write_alignment_file(
            trimmed_alignment, out_format, cutoff_trimmed_aln_fname, data_type
        )
        logging.info(
            'Wrote trimmed alignment {}\n'.format(cutoff_trimmed_aln_fname)
        )
    logging.info('Plotting distance distributions and cutoffs')
    plotting_wrapper(
        distances, window_size, method, criterion, cutoffs, manual_cutoffs
    )


def get_stride(window_size, overlap):
    """Convert window overlap to stride for use in slicing."""
    return window_size - overlap


def check_cutoff_value(criterion, cutoff_value):
    """Validate single cutoff value from config file."""
    if criterion == 'lognorm':
        if cutoff_value > 0 and cutoff_value < 1:
            pass
        else:
            exit(
                'Invalid lognorm quantile cutoff value "{}". The cutoffs must be numbers greater than 0 and less than 1.'.format(
                    cutoff_value
                )
            )
    elif criterion == 'mean':
        if cutoff_value > 1:
            pass
        elif cutoff_value > 0 and cutoff_value < 1:
            print(
                'WARNING: cutoff value "{}" is less than 1 mean. Did you intend to specify "lognorm" as criterion?'.format(
                    cutoff_value
                )
            )
        elif cutoff_value < 0:
            exit(
                'Invalid mean cutoff value "{}". Cutoffs must be greater than 0.'
            )


def check_cutoffs(criterion, cutoffs):
    """Wrapper to validate multiple cutoffs from config."""
    for cutoff in cutoffs:
        try:
            cutoff_value = float(cutoff)
            check_cutoff_value(criterion, cutoff_value)
        except ValueError as ex:
            exit(
                'Invalid cutoff value, cannot convert "{}" to number.'.format(
                    cutoff
                )
            )


def check_manual_cutoffs(criterion, manual_cutoffs):
    """Wrapper to validate manual cutoffs from config."""
    if manual_cutoffs:
        cutoffs = []
        for cutoff_tuple in manual_cutoffs:
            try:
                taxon, cutoff = cutoff_tuple
                cutoffs.append(cutoff)
            except ValueError as ex:
                exit(
                    'You set manual cutoffs that cannot be parsed. Make sure your configuration file is formatted correctly.'
                )
        check_cutoffs(criterion, cutoffs)


def get_validated_input(parsed_config):
    """Validate input from configuration file.

    Given parameters in configuration file either exit or return dict of vetted input.
    """
    valid_input_dict = {}
    # input config section
    try:
        alignment_name = parsed_config.get('input', 'input_file_name')
        with open(alignment_name) as f:
            pass
        valid_input_dict['alignment_name'] = alignment_name
    except IOError as ex:
        exit(
            'Sorry, could not read input alignment file "{}": {}'.format(
                alignment_name, ex.strerror
            )
        )
    file_format = parsed_config.get('input', 'input_format')
    if (
        file_format == 'fasta'
        or file_format == 'phylip'
        or file_format == 'phylip-int'
        or file_format == 'nexus'
        or file_format == 'nexus-int'
    ):
        valid_input_dict['file_format'] = file_format
    else:
        exit(
            'Invalid input file format: "{}". Choose from: fasta, phylip, phylip-int, nexus, or nexus-int.'.format(
                file_format
            )
        )
    data_type = parsed_config.get('input', 'data_type')
    if data_type == 'aa' or data_type == 'nt':
        valid_input_dict['data_type'] = data_type
    else:
        exit(
            'Invalid data type: "{}". Choose from: aa or nt.'.format(data_type)
        )
    try:
        distances_json = parsed_config.get('input', 'distances_object_file')
        if distances_json:
            with open(distances_json) as f:
                valid_input_dict['distances_json'] = distances_json
                pass
        else:
            valid_input_dict['distances_json'] = None
    except IOError as ex:
        exit(
            'Sorry, could not read distances file "{}": {}'.format(
                distances_json, ex.strerror
            )
        )
    try:
        tree_file = parsed_config.get('input', 'guide_tree')
        if tree_file:
            with open(tree_file) as f:
                pass
            valid_input_dict['tree_file'] = tree_file
        else:
            valid_input_dict['tree_file'] = None
    except IOError as ex:
        exit(
            'Sorry, could not read tree file "{}": {}'.format(
                tree_file, ex.strerror
            )
        )
    method = parsed_config.get('analysis', 'distance_method')
    if method == 'uncorrected' or method == 'jc':
        valid_input_dict['method'] = method
    else:
        exit(
            'Invalid distance method: "{}". Choose between "uncorrected" and "jc".'.format(
                method
            )
        )
    # analysis config section
    criterion = parsed_config.get('analysis', 'criterion')
    if criterion == 'lognorm' or criterion == 'mean':
        valid_input_dict['criterion'] = criterion
    else:
        exit(
            'Invalid criterion "{}". Choose between "lognorm" and "mean".'.format(
                criterion
            )
        )
    try:
        window_size = parsed_config.getint('analysis', 'window_size')
        overlap = parsed_config.getint('analysis', 'overlap')
        if window_size < 0:
            exit(
                'Invalid window size "{}". Window size cannot be smaller than 0.'.format(
                    window_size
                )
            )
        else:
            valid_input_dict['window_size'] = window_size
        if overlap >= 0 and overlap < window_size:
            valid_input_dict['overlap'] = overlap
        else:
            exit(
                'Invalid overlap "{}" for window size "{}". Overlap has to be integer smaller than window size.'.format(
                    overlap, window_size
                )
            )
        fraction = parsed_config.getfloat('analysis', 'fraction')
        if fraction >= 0 and fraction <= 1:
            valid_input_dict['fraction'] = fraction
        else:
            exit(
                'Invalid taxon fraction value "{}". Fraction must be between 0 and 1.'.format(
                    fraction
                )
            )
        cores = parsed_config.getint('analysis', 'cores')
        available_cpus = mp.cpu_count()
        if cores > available_cpus:
            exit(
                'You specified more ({}) compute cores than are available ({}). Exiting.'.format(
                    cores, available_cpus
                )
            )
        else:
            valid_input_dict['cores'] = cores
    except ValueError as ex:
        exit(
            'Invalid number input somewhere in your analysis configuration: {}.'.format(
                ex
            )
        )
    cutoffs = parsed_config.get('analysis', 'cutoffs').split(',')
    check_cutoffs(criterion, cutoffs)
    valid_input_dict['cutoffs'] = cutoffs
    manual_cutoffs = parsed_config.get('analysis', 'manual_cutoffs')
    if not manual_cutoffs:
        manual_cutoffs = None
    else:
        manual_cutoffs = [
            tuple(taxon_cutoff.split(','))
            for taxon_cutoff in parsed_config.get(
                'analysis', 'manual_cutoffs'
            ).split(';')
        ]
    check_manual_cutoffs(criterion, manual_cutoffs)
    valid_input_dict['manual_cutoffs'] = manual_cutoffs
    # output config section
    output_format = parsed_config.get('output', 'output_format')
    if (
        output_format == 'fasta'
        or output_format == 'phylip'
        or output_format == 'phylip-int'
        or output_format == 'nexus'
        or output_format == 'nexus-int'
    ):
        valid_input_dict['output_format'] = output_format
    else:
        exit(
            'Invalid output file format: "{}". Choose from: fasta, phylip, phylip-int, nexus, or nexus-int.'.format(
                output_format
            )
        )
    try:
        output_file_aln = parsed_config.get('output', 'output_file_aln')
        valid_input_dict['output_file_aln'] = output_file_aln
        report = parsed_config.get('output', 'report')
        valid_input_dict['report'] = report
        log_file_name = parsed_config.get('output', 'log')
        valid_input_dict['log_file_name'] = log_file_name
    except IOError as ex:
        exit('Invalid output: {}'.format(ex.strerror))
    return valid_input_dict


def main():
    start = time.time()
    try:
        script, config_file_name = argv
    except ValueError as ex:
        exit('You need to provide a configuration file.\nUsage: spruceup.py <my_config_file>')
    conf = read_config(config_file_name)
    valid_input_dict = get_validated_input(conf)
    level = logging.INFO
    log_format = '%(asctime)s - %(message)s'
    handlers = [
        logging.FileHandler(valid_input_dict['log_file_name'], 'w'),
        logging.StreamHandler(),
    ]
    logging.basicConfig(level=level, format=log_format, handlers=handlers)
    # check for existing distances file
    if valid_input_dict['distances_json']:
        logging.info(
            'Parsing alignment {} ...\n'.format(
                valid_input_dict['alignment_name']
            )
        )
        alignment = aln_parsing.parse_alignment(
            valid_input_dict['alignment_name'], valid_input_dict['file_format']
        )
        mean_taxon_distances = read_distances_dict(
            valid_input_dict['distances_json']
        )
    else:
        alignment, mean_taxon_distances = analyze(
            valid_input_dict['alignment_name'],
            valid_input_dict['tree_file'],
            valid_input_dict['file_format'],
            valid_input_dict['data_type'],
            valid_input_dict['window_size'],
            valid_input_dict['overlap'],
            valid_input_dict['cores'],
            valid_input_dict['method'],
            valid_input_dict['fraction'],
        )
    print_mem()
    output_loop(
        alignment,
        mean_taxon_distances,
        valid_input_dict['window_size'],
        valid_input_dict['method'],
        valid_input_dict['criterion'],
        valid_input_dict['cutoffs'],
        valid_input_dict['manual_cutoffs'],
        valid_input_dict['report'],
        valid_input_dict['output_format'],
        valid_input_dict['output_file_aln'],
        valid_input_dict['data_type'],
    )
    logging.info('Finished in {:.2f} seconds'.format(time.time() - start))


if __name__ == '__main__':
    
    main()
