#! /usr/bin/env python3

from statistics import mean
from math import log

import unittest

from spruceup import spruceup

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


class TestGetDistances(unittest.TestCase):
    def test_get_uncorrected_no_nt_distances_returns_both_aln_name_dists(self):
        result = len(
            spruceup.get_distances(
                alns_tpls[0],
                tree_dists,
                methods[0],
                fractions[0],
                data_types[1],
            )
        )
        self.assertEqual(result, 2)

    def test_get_uncorrected_half_nt_distances_returns_both_aln_name_dists(
        self
    ):
        result = len(
            spruceup.get_distances(
                alns_tpls[0],
                tree_dists,
                methods[0],
                fractions[1],
                data_types[1],
            )
        )
        self.assertEqual(result, 2)

    def test_get_uncorrected_all_aa_distances_returns_both_aln_name_dists(self):
        result = len(
            spruceup.get_distances(
                alns_tpls[0],
                tree_dists,
                methods[1],
                fractions[2],
                data_types[0],
            )
        )
        self.assertEqual(result, 2)

    def test_get_jc_all_nt_distances_returns_both_aln_name_dists(self):
        result = len(
            spruceup.get_distances(
                alns_tpls[0],
                tree_dists,
                methods[1],
                fractions[2],
                data_types[1],
            )
        )
        self.assertEqual(result, 2)

    def test_get_uncorrected_all_nt_distances_returns_both_aln_name_dists(self):
        result = len(
            spruceup.get_distances(
                alns_tpls[0],
                tree_dists,
                methods[0],
                fractions[2],
                data_types[1],
            )
        )
        self.assertEqual(result, 2)

    def test_get_uncorrected_nt_distances_returns_zero_dist_on_ident_seqs(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[0], tree_dists, methods[0], fractions[2], data_types[1]
        )
        dists_list = []
        for dist_tpl in distances:
            sp1, sp2, dist = dist_tpl
            dists_list.append(dist)
        avg_dist = spruceup.get_list_mean(dists_list)
        self.assertEqual(avg_dist, 0)

    def test_get_uncorrected_nt_distances_returns_non_zero_dist_on_different_seqs1(
        self
    ):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[1], tree_dists, methods[0], fractions[2], data_types[1]
        )
        dists_list = []
        for dist_tpl in distances:
            sp1, sp2, dist = dist_tpl
            dists_list.append(dist)
        avg_dist = spruceup.get_list_mean(dists_list)
        self.assertFalse(avg_dist == 0)

    def test_get_uncorrected_nt_distances_returns_non_zero_dist_on_different_seqs2(
        self
    ):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        dists_list = []
        for dist_tpl in distances:
            sp1, sp2, dist = dist_tpl
            dists_list.append(dist)
        avg_dist = spruceup.get_list_mean(dists_list)
        self.assertFalse(avg_dist == 0)


class TestGetDistAndTaxa(unittest.TestCase):
    def test_find_all_dists_on_simple_aln(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        self.assertEqual(len(distances), len(dists))

    def test_find_all_taxa_lists_on_simple_aln(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        self.assertTrue(set(spp).issubset(taxa_rows))


class TestDistMatrix(unittest.TestCase):
    def test_find_all_taxa_lists_on_simple_aln(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        self.assertTrue(set(spp).issubset(taxa_rows))

    def test_find_all_dists_on_simple_aln(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        self.assertEqual(len(distances), len(dists))

    def test_dist_matrix_dimensions_on_simple_aln(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        aln_name, aln_dict = alns_tpls[2]
        no_taxa = len(aln_dict.keys())
        matrix = spruceup.get_dist_matrix(dists, no_taxa)
        self.assertEqual(len(matrix), no_taxa)
        self.assertEqual(len(matrix[0]), no_taxa)

    def test_taxon_map_returns_all_taxa(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[2], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        aln_name, aln_dict = alns_tpls[2]
        taxa = aln_dict.keys()
        taxon_map = spruceup.get_taxon_map(taxa_rows)
        self.assertTrue(set(taxa).issubset(taxa_rows))

    def test_mean_distances_returns_all_zeros_if_seqs_identical(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[0], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        aln_name, aln_dict = alns_tpls[0]
        no_taxa = len(aln_dict.keys())
        matrix = spruceup.get_dist_matrix(dists, no_taxa)
        taxon_map = spruceup.get_taxon_map(taxa_rows)
        mean_distances = spruceup.get_mean_distances(matrix, taxon_map)
        self.assertTrue(
            all(mean_dist == 0 for mean_dist in mean_distances.values())
        )

    def test_mean_distances_returns_non_zeros_if_seqs_different(self):
        aln_name, distances = spruceup.get_distances(
            alns_tpls[1], tree_dists, methods[0], fractions[2], data_types[1]
        )
        taxa_rows, dists = spruceup.get_dist_and_taxa_lists(distances)
        aln_name, aln_dict = alns_tpls[1]
        no_taxa = len(aln_dict.keys())
        matrix = spruceup.get_dist_matrix(dists, no_taxa)
        taxon_map = spruceup.get_taxon_map(taxa_rows)
        mean_distances = spruceup.get_mean_distances(matrix, taxon_map)
        self.assertFalse(
            all(mean_dist == 0 for mean_dist in mean_distances.values())
        )


class TestRanAndScaledDists(unittest.TestCase):
    def test_jc_nt_returns_zero_if_dist_is_zero(self):
        self.assertEqual(
            spruceup.jc_correction((10, 0), data_types[1]), (10, 0)
        )

    def test_jc_returns_nan_if_dist_is_nan(self):
        self.assertEqual(
            spruceup.jc_correction((10, 'NaN'), data_types[1]), (10, 'NaN')
        )

    def test_jc_returns_zero_if_nt_seqs_identical(self):
        self.assertEqual(
            spruceup.jc_correction(
                spruceup.p_distance(simple_seq, simple_seq), data_types[1]
            ),
            (10, 0),
        )

    def test_jc_returns_zero_if_aa_seqs_identical(self):
        self.assertEqual(
            spruceup.jc_correction(
                spruceup.p_distance(simple_aa_seq, simple_aa_seq), data_types[0]
            ),
            (10, 0),
        )

    def test_jc_returns_max_nt_dist_if_10_out_of_10_nt_diff(self):
        self.assertEqual(
            spruceup.jc_correction((10, 10), data_types[1]),
            (10, (3 / 4 * log(1 - 4 / 3 * -10))),
        )

    def test_jc_returns_max_nt_dist_if_max_diff_nt(self):
        dist = spruceup.p_distance(simple_seq, different_seq)
        self.assertEqual(
            spruceup.jc_correction(dist, data_types[1]),
            (dist[0], (3 / 4 * log(1 - 4 / 3 * -dist[1]))),
        )

    def test_jc_returns_max_aa_dist_if_10_out_of_10_nt_diff(self):
        self.assertEqual(
            spruceup.jc_correction((10, 10), data_types[0]),
            (10, (19 / 20 * log(1 - 20 / 19 * -10))),
        )

    def test_jc_returns_max_aa_dist_if_max_diff_nt(self):
        dist = spruceup.p_distance(simple_aa_seq, different_aa_seq)
        self.assertEqual(
            spruceup.jc_correction(dist, data_types[0]),
            (dist[0], (19 / 20 * log(1 - 20 / 19 * -dist[1]))),
        )

    def test_scaled_returns_nan_if_one_seq_empty(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(simple_seq, empty_seq)
            ),
            'NaN',
        )

    def test_scaled_returns_nan_if_both_seqs_empty(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(empty_seq, empty_seq)
            ),
            'NaN',
        )

    def test_scaled_returns_1_if_max_seq_difference(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(simple_seq, different_seq)
            ),
            1,
        )

    def test_scaled_returns_01_if_10percent_seq_different(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(simple_seq, partly_gappy)
            ),
            0.1,
        )

    def test_scaled_returns_0_if_seqs_identical1(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(partly_gappy, partly_gappy)
            ),
            0.0,
        )

    def test_scaled_returns_0_if_seqs_identical2(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(simple_seq, simple_seq)
            ),
            0.0,
        )

    def test_scaled_returns_0_if_seqs_identical3(self):
        self.assertEqual(
            spruceup.get_scaled_distance(
                spruceup.p_distance(simple_aa_seq, simple_aa_seq)
            ),
            0.0,
        )

    def test_aa_distance_if_both_sequences_max_different(self):
        self.assertEqual(
            spruceup.p_distance(simple_aa_seq, different_aa_seq), (10, 10)
        )

    def test_nt_distance_if_both_sequences_max_different(self):
        self.assertEqual(
            spruceup.p_distance(simple_seq, different_seq), (10, 10)
        )

    def test_returns_efflen0_and_nan_if_empty_seq_first(self):
        self.assertEqual(spruceup.p_distance(empty_seq, simple_seq), (0, 'NaN'))

    def test_returns_efflen10_and_nan_if_empty_seq_second(self):
        self.assertEqual(
            spruceup.p_distance(simple_seq, empty_seq), (10, 'NaN')
        )

    def test_returns_efflen7_and_1differnce_if_partly_gappy_seq_first(self):
        self.assertEqual(spruceup.p_distance(partly_gappy, simple_seq), (7, 1))

    def test_returns_efflen10_and_1differnce_if_partly_gappy_seq_second(self):
        self.assertEqual(spruceup.p_distance(simple_seq, partly_gappy), (10, 1))

    def test_raises_error_if_seqs_unaligned(self):
        self.assertRaises(
            ValueError, spruceup.p_distance, simple_seq, short_seq
        )


if __name__ == '__main__':
    unittest.main()
