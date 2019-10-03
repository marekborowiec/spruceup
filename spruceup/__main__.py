import argparse
import logging
import os
import time

import psutil

from spruceup import aln_parsing, aln_writing, spruceup
from spruceup.__init__ import __version__

parser = argparse.ArgumentParser(usage='spruceup.py <config_file_name>')

parser.add_argument(
    'config_file_name',
    type = str,
    help = 'Configuration file with spruceup run parameters.'
    )
parser.add_argument(
    '-v',
    '--version', 
    action='version',
    version='%(prog)s {version}'.format(version=__version__),
    help='show the version number and exit'
    )


def print_mem():
    """Print memory usage using psutil module."""
    process = psutil.Process(os.getpid())
    mem = process.memory_info().rss / 1_000_000
    logging.info('Used {0:.2f} MB memory\n'.format(mem))


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
    no_missing_ambiguous_dict = spruceup.replace_missing_in_dict(aln_dict, data_type)
    if tree_file_name == None:
        logging.info('Did not find guide tree, continuing ...\n')
        tree_dists = None
    else:
        logging.info('Reading in guide tree {} ...\n'.format(tree_file_name))
        tree_dists = spruceup.get_tree_dist_dict(tree_file_name)
    windows = spruceup.get_windows(no_missing_ambiguous_dict, window_size, overlap)
    all_distances = spruceup.distances_wrapper(
        windows, tree_dists, cores, data_type, method=method, fraction=fraction
    )
    taxa_distances = spruceup.dist_taxa_wrapper(all_distances)
    mean_aln_distances = spruceup.mean_distances_wrapper(taxa_distances)
    mean_taxon_distances = spruceup.dists_per_taxon(mean_aln_distances)
    spruceup.write_distances_dict(mean_taxon_distances, method, window_size, overlap)
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
    alignment_sites = spruceup.get_alignment_size(untrimmed_alignment)
    cutoff_floats = [float(cutoff_string) for cutoff_string in cutoffs]
    for cutoff in cutoff_floats:
        logging.info(
            'Finding outliers for {} {}s cutoff ...'.format(cutoff, criterion)
        )
        outliers = spruceup.get_outliers_wrapper(
            distances, window_size, method, criterion, cutoff, manual_cutoffs
        )
        sites_removed, trimmed_alignment = spruceup.remove_outliers(
            untrimmed_alignment, outliers
        )
        logging.info('Sites removed: {}'.format(sites_removed))
        percent_removed = (
            spruceup.get_removed_fraction(alignment_sites, sites_removed) * 100
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
        spruceup.write_report(
            spruceup.print_report(outliers, criterion, cutoff, manual_cutoffs),
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
    spruceup.plotting_wrapper(
        distances, window_size, method, criterion, cutoffs, manual_cutoffs
    )


def main():
    start = time.time()
    args = parser.parse_args()
    config_file_name = args.config_file_name
    conf = spruceup.read_config(config_file_name)
    valid_input_dict = spruceup.get_validated_input(conf)
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
        mean_taxon_distances = spruceup.read_distances_dict(
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