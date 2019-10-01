#! /usr/bin/env python3

from statistics import mean
from math import log

from nose.tools import assert_true, assert_false, assert_equal, assert_raises

from spruceup import (
    jc_correction,
    p_distance,
    get_scaled_distance,
    get_distances,
    get_taxon_map,
    get_list_mean,
    get_dist_and_taxa_lists,
    get_dist_matrix,
    get_taxon_map,
    get_mean_distances,
)

# some short mock sequences
empty_seq = '-' * 10
partly_gappy = '---CAGTACC'
simple_seq = 'TGTCAGTACG'
different_seq = 'GCAATAATAA'
short_seq = 'CAGTCG'
simple_aa_seq = 'LNIWVNLAVC'
different_aa_seq = 'DKVCTIVSTQ'

# different dicts represent different parsed alignments
aln_dict1 = {
    'sp1': simple_seq,
    'sp2': simple_seq,
    'sp3': simple_seq,
    'sp4': simple_seq,
    'sp5': simple_seq,
}
aln_dict2 = {
    'sp1': simple_seq,
    'sp2': different_seq,
    'sp3': empty_seq,
    'sp4': empty_seq,
    'sp5': empty_seq,
}
aln_dict3 = {
    'sp1': simple_seq,
    'sp2': different_seq,
    'sp3': simple_seq,
    'sp4': partly_gappy,
    'sp5': empty_seq,
}

aln_ident = ('aln_ident', aln_dict1)
aln_gappy = ('aln_gappy', aln_dict2)
aln_norma = ('aln_norma', aln_dict3)
alns_tpls = [aln_ident, aln_gappy, aln_norma]

spp = ['sp1', 'sp2', 'sp3', 'sp4', 'sp5']
zero_dist_list = [(sp1, sp2, 0) for sp1 in spp for sp2 in spp]

methods = ['uncorrected', 'jc']

fractions = [0, 0.5, 1]

data_types = ['aa', 'nt']

tree_dists = {
    ('sp1', 'sp1'): 0.0,
    ('sp1', 'sp2'): 0.1,
    ('sp1', 'sp3'): 3,
    ('sp1', 'sp4'): 0.5,
    ('sp1', 'sp5'): 0.0,
    ('sp2', 'sp1'): 0.1,
    ('sp2', 'sp2'): 0.0,
    ('sp2', 'sp3'): 0.5,
    ('sp2', 'sp4'): 0.1,
    ('sp2', 'sp5'): 0.1,
    ('sp3', 'sp1'): 3,
    ('sp3', 'sp2'): 2.5,
    ('sp3', 'sp3'): 0.0,
    ('sp3', 'sp4'): 0.9,
    ('sp3', 'sp5'): 1,
    ('sp4', 'sp1'): 0.5,
    ('sp4', 'sp2'): 0.1,
    ('sp4', 'sp3'): 0.9,
    ('sp4', 'sp4'): 0.0,
    ('sp4', 'sp5'): 0.7,
    ('sp5', 'sp1'): 0.0,
    ('sp5', 'sp2'): 0.1,
    ('sp5', 'sp3'): 1,
    ('sp5', 'sp4'): 0.7,
    ('sp5', 'sp5'): 0.0,
}


def test_get_uncorrected_no_nt_distances_returns_both_aln_name_dists():
    assert_equal(
        len(
            get_distances(
                alns_tpls[0],
                tree_dists,
                methods[0],
                fractions[0],
                data_types[1],
            )
        ),
        2,
    )


def test_get_uncorrected_half_nt_distances_returns_both_aln_name_dists():
    assert_equal(
        len(
            get_distances(
                alns_tpls[0],
                tree_dists,
                methods[0],
                fractions[1],
                data_types[1],
            )
        ),
        2,
    )


def test_get_uncorrected_all_aa_distances_returns_both_aln_name_dists():
    assert_equal(
        len(
            get_distances(
                alns_tpls[0],
                tree_dists,
                methods[1],
                fractions[2],
                data_types[0],
            )
        ),
        2,
    )


def test_get_jc_all_nt_distances_returns_both_aln_name_dists():
    assert_equal(
        len(
            get_distances(
                alns_tpls[0],
                tree_dists,
                methods[1],
                fractions[2],
                data_types[1],
            )
        ),
        2,
    )


def test_get_uncorrected_all_nt_distances_returns_both_aln_name_dists():
    assert_equal(
        len(
            get_distances(
                alns_tpls[0],
                tree_dists,
                methods[0],
                fractions[2],
                data_types[1],
            )
        ),
        2,
    )


def test_get_uncorrected_nt_distances_returns_zero_dist_on_ident_seqs():
    aln_name, distances = get_distances(
        alns_tpls[0], tree_dists, methods[0], fractions[2], data_types[1]
    )
    dists_list = []
    for dist_tpl in distances:
        sp1, sp2, dist = dist_tpl
        dists_list.append(dist)
    avg_dist = get_list_mean(dists_list)
    assert_equal(avg_dist, 0)


def test_get_uncorrected_nt_distances_returns_non_zero_dist_on_different_seqs1():
    aln_name, distances = get_distances(
        alns_tpls[1], tree_dists, methods[0], fractions[2], data_types[1]
    )
    dists_list = []
    for dist_tpl in distances:
        sp1, sp2, dist = dist_tpl
        dists_list.append(dist)
    print(dists_list)
    avg_dist = get_list_mean(dists_list)
    assert_false(avg_dist == 0)


def test_get_uncorrected_nt_distances_returns_non_zero_dist_on_different_seqs2():
    aln_name, distances = get_distances(
        alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
    )
    dists_list = []
    for dist_tpl in distances:
        sp1, sp2, dist = dist_tpl
        dists_list.append(dist)
    print(dists_list)
    avg_dist = get_list_mean(dists_list)
    assert_false(avg_dist == 0)


def test_find_all_taxa_lists_on_simple_aln():
    aln_name, distances = get_distances(
        alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
    )
    taxa_rows, dists = get_dist_and_taxa_lists(distances)
    assert_true(set(spp).issubset(taxa_rows))


def test_find_all_dists_on_simple_aln():
    aln_name, distances = get_distances(
        alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
    )
    taxa_rows, dists = get_dist_and_taxa_lists(distances)
    assert_equal(len(distances), len(dists))


def test_dist_matrix_dimensions_on_simple_aln():
    aln_name, distances = get_distances(
        alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
    )
    taxa_rows, dists = get_dist_and_taxa_lists(distances)
    aln_name, aln_dict = alns_tpls[2]
    no_taxa = len(aln_dict.keys())
    matrix = get_dist_matrix(dists, no_taxa)
    assert_equal(len(matrix), no_taxa)
    assert_equal(len(matrix[0]), no_taxa)


def test_taxon_map_returns_all_taxa():
    aln_name, distances = get_distances(
        alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
    )
    taxa_rows, dists = get_dist_and_taxa_lists(distances)
    aln_name, aln_dict = alns_tpls[2]
    taxa = aln_dict.keys()
    taxon_map = get_taxon_map(taxa_rows)
    assert_true(set(taxa).issubset(taxa_rows))


def test_mean_distances_returns_all_zeros_if_seqs_identical():
    aln_name, distances = get_distances(
        alns_tpls[0], tree_dists, methods[0], fractions[2], data_types[1]
    )
    taxa_rows, dists = get_dist_and_taxa_lists(distances)
    aln_name, aln_dict = alns_tpls[0]
    no_taxa = len(aln_dict.keys())
    matrix = get_dist_matrix(dists, no_taxa)
    taxon_map = get_taxon_map(taxa_rows)
    mean_distances = get_mean_distances(matrix, taxon_map)
    assert_true(all(mean_dist == 0 for mean_dist in mean_distances.values()))


def test_mean_distances_returns_non_zeros_if_seqs_different():
    aln_name, distances = get_distances(
        alns_tpls[1], tree_dists, methods[0], fractions[2], data_types[1]
    )
    taxa_rows, dists = get_dist_and_taxa_lists(distances)
    aln_name, aln_dict = alns_tpls[1]
    no_taxa = len(aln_dict.keys())
    matrix = get_dist_matrix(dists, no_taxa)
    taxon_map = get_taxon_map(taxa_rows)
    mean_distances = get_mean_distances(matrix, taxon_map)
    assert_false(all(mean_dist == 0 for mean_dist in mean_distances.values()))


def test_jc_nt_returns_zero_if_dist_is_zero():
    assert_equal(jc_correction((10, 0), data_types[1]), (10, 0))


def test_jc_returns_nan_if_dist_is_nan():
    assert_equal(jc_correction((10, 'NaN'), data_types[1]), (10, 'NaN'))


def test_jc_returns_zero_if_nt_seqs_identical():
    assert_equal(
        jc_correction(p_distance(simple_seq, simple_seq), data_types[1]),
        (10, 0),
    )


def test_jc_returns_zero_if_aa_seqs_identical():
    assert_equal(
        jc_correction(p_distance(simple_aa_seq, simple_aa_seq), data_types[0]),
        (10, 0),
    )


def test_jc_returns_max_nt_dist_if_10_out_of_10_nt_diff():
    assert_equal(
        jc_correction((10, 10), data_types[1]),
        (10, (3 / 4 * log(1 - 4 / 3 * -10))),
    )


def test_jc_returns_max_nt_dist_if_max_diff_nt():
    dist = p_distance(simple_seq, different_seq)
    assert_equal(
        jc_correction(dist, data_types[1]),
        (dist[0], (3 / 4 * log(1 - 4 / 3 * -dist[1]))),
    )


def test_jc_returns_max_aa_dist_if_10_out_of_10_nt_diff():
    assert_equal(
        jc_correction((10, 10), data_types[0]),
        (10, (19 / 20 * log(1 - 20 / 19 * -10))),
    )


def test_jc_returns_max_aa_dist_if_max_diff_nt():
    dist = p_distance(simple_aa_seq, different_aa_seq)
    assert_equal(
        jc_correction(dist, data_types[0]),
        (dist[0], (19 / 20 * log(1 - 20 / 19 * -dist[1]))),
    )


def test_scaled_returns_nan_if_one_seq_empty():
    assert_equal(get_scaled_distance(p_distance(simple_seq, empty_seq)), 'NaN')


def test_scaled_returns_nan_if_both_seqs_empty():
    assert_equal(get_scaled_distance(p_distance(empty_seq, empty_seq)), 'NaN')


def test_scaled_returns_1_if_max_seq_difference():
    assert_equal(get_scaled_distance(p_distance(simple_seq, different_seq)), 1)


def test_scaled_returns_01_if_10percent_seq_different():
    assert_equal(get_scaled_distance(p_distance(simple_seq, partly_gappy)), 0.1)


def test_scaled_returns_0_if_seqs_identical1():
    assert_equal(
        get_scaled_distance(p_distance(partly_gappy, partly_gappy)), 0.0
    )


def test_scaled_returns_0_if_seqs_identical2():
    assert_equal(get_scaled_distance(p_distance(simple_seq, simple_seq)), 0.0)


def test_scaled_returns_0_if_seqs_identical3():
    assert_equal(
        get_scaled_distance(p_distance(simple_aa_seq, simple_aa_seq)), 0.0
    )


def test_aa_distance_if_both_sequences_max_different():
    assert_equal(p_distance(simple_aa_seq, different_aa_seq), (10, 10))


def test_nt_distance_if_both_sequences_max_different():
    assert_equal(p_distance(simple_seq, different_seq), (10, 10))


def test_returns_efflen0_and_nan_if_empty_seq_first():
    assert_equal(p_distance(empty_seq, simple_seq), (0, 'NaN'))


def test_returns_efflen10_and_nan_if_empty_seq_second():
    assert_equal(p_distance(simple_seq, empty_seq), (10, 'NaN'))


def test_returns_efflen7_and_1differnce_if_partly_gappy_seq_first():
    assert_equal(p_distance(partly_gappy, simple_seq), (7, 1))


def test_returns_efflen10_and_1differnce_if_partly_gappy_seq_second():
    assert_equal(p_distance(simple_seq, partly_gappy), (10, 1))


def test_raises_error_if_seqs_unaligned():
    assert_raises(ValueError, p_distance, simple_seq, short_seq)

if __name__ == '__main__':
    test_get_uncorrected_no_nt_distances_returns_both_aln_name_dists()
    test_get_uncorrected_half_nt_distances_returns_both_aln_name_dists()
    test_get_uncorrected_all_aa_distances_returns_both_aln_name_dists()
    test_get_jc_all_nt_distances_returns_both_aln_name_dists()
    test_get_uncorrected_all_nt_distances_returns_both_aln_name_dists()
    test_get_uncorrected_nt_distances_returns_zero_dist_on_ident_seqs()
    test_get_uncorrected_nt_distances_returns_non_zero_dist_on_different_seqs1()
    test_get_uncorrected_nt_distances_returns_non_zero_dist_on_different_seqs2()
    test_find_all_taxa_lists_on_simple_aln()
    test_find_all_dists_on_simple_aln()
    test_dist_matrix_dimensions_on_simple_aln()
    test_taxon_map_returns_all_taxa()
    test_mean_distances_returns_all_zeros_if_seqs_identical()
    test_mean_distances_returns_non_zeros_if_seqs_different()
    test_jc_nt_returns_zero_if_dist_is_zero()
    test_jc_returns_nan_if_dist_is_nan()
    test_jc_returns_zero_if_nt_seqs_identical()
    test_jc_returns_zero_if_aa_seqs_identical()
    test_jc_returns_max_nt_dist_if_10_out_of_10_nt_diff()
    test_jc_returns_max_nt_dist_if_max_diff_nt()
    test_jc_returns_max_aa_dist_if_10_out_of_10_nt_diff()
    test_jc_returns_max_aa_dist_if_max_diff_nt()
    test_scaled_returns_nan_if_one_seq_empty()
    test_scaled_returns_nan_if_both_seqs_empty()
    test_scaled_returns_1_if_max_seq_difference()
    test_scaled_returns_01_if_10percent_seq_different()
    test_scaled_returns_0_if_seqs_identical1()
    test_scaled_returns_0_if_seqs_identical2()
    test_scaled_returns_0_if_seqs_identical3()
    test_aa_distance_if_both_sequences_max_different()
    test_nt_distance_if_both_sequences_max_different()
    test_returns_efflen0_and_nan_if_empty_seq_first()
    test_returns_efflen10_and_nan_if_empty_seq_second()
    test_returns_efflen7_and_1differnce_if_partly_gappy_seq_first()
    test_returns_efflen10_and_1differnce_if_partly_gappy_seq_second()
    test_raises_error_if_seqs_unaligned()