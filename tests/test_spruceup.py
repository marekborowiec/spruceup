from math import log
import random

from nose.tools import assert_equal, assert_raises

def p_distance(seq1, seq2):
    """Calculate Hamming/p-distance for two sequences.

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
            if el1 is not '-' and el2 is not '-' and el1 is not '?' and el2 is not '?'
        )
    else:
        p_distance = 'NaN'
    return (eff_len1, p_distance)


def get_scaled_distance(distance_tpl):
    """Scale distances from 0 to 1.

    Given tuple of (seq length, Hamming distance)
    return proportion of different sites
    unless one of the sequences was missing."""
    eff_seq_len, distance = distance_tpl
    if distance is not 'NaN':
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
    return tuple of (seq length, corrected distance)."""
    eff_seq_len, p_distance = distance_tpl
    if p_distance is not 'NaN':
        if data_type == 'nt':
            jc_corrected = 3 / 4 * log(1 - 4 / 3 * -p_distance)
        elif data_type == 'aa':
            jc_corrected = 19 / 20 * log(1 - 20 / 19 * -p_distance)
    else:
            jc_corrected = 'NaN'
    return (eff_seq_len, jc_corrected)

def get_distances(aln_tuple, method, fraction, data_type):
    """Calculate uncorrected or JC-corrected distances for alignment.

    Given tuple (alignment name, alignment distances dict),
    return tuple of (alignment name, list of pairwise distances).
    Args:
    method (str) -- 'uncorrected' or 'jc'
    fraction (int) -- 0 to 1
    data_type (str) -- 'nt' or 'aa'
    """

    aln_name, aln_dict = aln_tuple
    seqs_to_compare_to = random.sample(
        aln_dict.items(), int(len(aln_dict.items()) * fraction)
    )

    if method == 'uncorrected':
        distances = [
            (sp1, sp2, get_scaled_distance(p_distance(seq1, seq2)))
            for sp2, seq2 in seqs_to_compare_to
            for sp1, seq1 in aln_dict.items()
            if sp1 != sp2
        ]
    elif method == 'jc':
        distances = [
            (sp1, sp2, get_scaled_distance(jc_correction(p_distance(seq1, seq2), data_type)))
            for sp2, seq2 in seqs_to_compare_to
            for sp1, seq1 in aln_dict.items()
            if sp1 != sp2
        ]
    return (aln_name, distances)


empty_seq = '-' * 10
partly_gappy = '---CAGTACC'
simple_seq = 'TGTCAGTACG'
different_seq = 'GCAATAATAA'
short_seq = 'CAGTCG'
simple_aa_seq = 'LNIWVNLAVC'
different_aa_seq = 'DKVCTIVSTQ'
aln_dict1 = {'sp1' : simple_seq,
             'sp2' : simple_seq,
             'sp3' : simple_seq,
             'sp4' : simple_seq,
             'sp5' : simple_seq}
aln_dict2 = {'sp1' : simple_seq,
             'sp2' : empty_seq,
             'sp3' : empty_seq,
             'sp4' : empty_seq,
             'sp5' : empty_seq}
aln_dict3 = {'sp1' : simple_seq,
             'sp2' : different_seq,
             'sp3' : simple_seq,
             'sp4' : partly_gappy,
             'sp5' : simple_seq}

aln_ident = ('aln_ident', aln_dict1)
aln_gappy = ('aln_gappy', aln_dict2)
aln_norma = ('aln_norma', aln_dict3)

spp = ['sp1', 'sp2', 'sp3', 'sp4', 'sp5']
zero_dist_list = [(sp1, sp2, 0) for sp1 in spp for sp2 in spp]

ident_dists = get_distances(aln_ident, 'jc', 1, 'nt')
gappy_dists = get_distances(aln_gappy, 'jc', 1, 'nt')
norma_dists = get_distances(aln_norma, 'uncorrected', 1, 'nt')

def get_dist_and_taxa_lists(distances):
    """Get aligned lists of taxa and distances.

    Given pairwise distances tuples of 
    (taxon1, taxon2, pairwise distance between the two)
    return tuple of (taxa, distances)."""
    taxa_rows = [sp2 for (sp1, sp2, dist) in distances]
    dists = [dist for (sp1, sp2, dist) in distances]
    return (taxa_rows, dists)

def dist_taxa_wrapper(dist_tuples):
    """Wrapper for getting aligned lists of taxa and distances.

    Given list of tuples of multiple alignments
    return tuple (alignment name : (taxa rows, distances list))."""
    for aln_name, distances in dist_tuples:
        return (aln_name, get_dist_and_taxa_lists(distances))


ident_l = dist_taxa_wrapper([ident_dists])
gappy_l = dist_taxa_wrapper([gappy_dists])
norma_l = dist_taxa_wrapper([norma_dists])

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
        #print(taxon_map)
        dist_matrix = get_dist_matrix(dists, len(taxon_map))
        print(dist_matrix)
        means = get_mean_distances(dist_matrix, taxon_map)
        #print(means)
        return (aln_name, means)


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
    list_mean = sum([i for i in lst if i is not 'NaN']) / float(len(lst))
    #print(lst)
    return round(list_mean, 5)

print(mean_distances_wrapper([gappy_l]))

"""
How should this work?
def test_ident_distances():
    assert_equal(get_distances(aln_ident, 'jc', 1, 'nt'),  ('aln_ident', [set(('sp1', 'sp2', 0), ('sp1', 'sp3', 0), 
                                                           ('sp1', 'sp4', 0), ('sp1', 'sp5', 0), ('sp2', 'sp1', 0), 
                                                           ('sp2', 'sp3', 0), ('sp2', 'sp4', 0), ('sp2', 'sp5', 0), 
                                                           ('sp3', 'sp1', 0), ('sp3', 'sp2', 0), ('sp3', 'sp4', 0), 
                                                           ('sp3', 'sp5', 0), ('sp4', 'sp1', 0), ('sp4', 'sp2', 0), 
                                                           ('sp4', 'sp3', 0), ('sp4', 'sp5', 0), ('sp5', 'sp1', 0), 
                                                           ('sp5', 'sp2', 0), ('sp5', 'sp3', 0), ('sp5', 'sp4', 0))]
    ))
"""
def test_jc_nt_zero():
    assert_equal(jc_correction((10, 0), 'nt'), (10, 0))

def test_jc_nan():
    assert_equal(jc_correction((10, 'NaN'), 'nt'), (10, 'NaN'))

def test_jc_nt_min():
    assert_equal(jc_correction(p_distance(simple_seq, simple_seq), 'nt'), (10, 0))

def test_jc_aa_min():
    assert_equal(jc_correction(p_distance(simple_aa_seq, simple_aa_seq), 'aa'), (10, 0))

def test_jc_nt_max():
    assert_equal(jc_correction((10, 10), 'nt'), (10, (3 / 4 * log(1 - 4 / 3 * -10))))

def test_scaled_nan():
    assert_equal(get_scaled_distance(p_distance(simple_seq, empty_seq)), 'NaN')


def test_scaled_max():
    assert_equal(get_scaled_distance(p_distance(simple_seq, different_seq)), 1)


def test_scaled_partial():
    assert_equal(get_scaled_distance(p_distance(simple_seq, partly_gappy)), 0.1)


def test_scaled_identical():
    assert_equal(get_scaled_distance(p_distance(partly_gappy, partly_gappy)), 0.0)

def test_aa_distance():
    assert_equal(p_distance(simple_aa_seq, different_aa_seq), (10, 10))


def test_nt_distance():
    assert_equal(p_distance(simple_seq, different_seq), (10, 10))


def test_gappy():
    assert_equal(p_distance(empty_seq, simple_seq), (0, 'NaN'))


def test_gappy_reverse():
    assert_equal(p_distance(simple_seq, empty_seq), (10, 'NaN'))


def test_partly_gappy():
    assert_equal(p_distance(partly_gappy, simple_seq), (7, 1))


def test_partly_gappy_reverse():
    assert_equal(p_distance(simple_seq, partly_gappy), (10, 1))

def test_unaligned():
    assert_raises(ValueError, p_distance, simple_seq, short_seq)