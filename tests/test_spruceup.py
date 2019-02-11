#! /usr/bin/env python3

# coding: utf-8
from math import log

from nose.tools import assert_equal, assert_raises

import aln_parsing
import aln_writing
from spruceup import jc_correction, p_distance, get_scaled_distance, get_distances, get_taxon_map

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
             'sp2' : different_seq,
             'sp3' : empty_seq,
             'sp4' : empty_seq,
             'sp5' : empty_seq}
aln_dict3 = {'sp1' : simple_seq,
             'sp2' : different_seq,
             'sp3' : simple_seq,
             'sp4' : partly_gappy,
             'sp5' : empty_seq}

aln_ident = ('aln_ident', aln_dict1)
aln_gappy = ('aln_gappy', aln_dict2)
aln_norma = ('aln_norma', aln_dict3)
alns_tpls = [aln_ident, aln_gappy, aln_norma]

spp = ['sp1', 'sp2', 'sp3', 'sp4', 'sp5']
zero_dist_list = [(sp1, sp2, 0) for sp1 in spp for sp2 in spp]

methods = ['uncorrected', 'jc']

fractions = [0, 0.5, 1]

data_types = ['aa', 'nt']

def test_get_taxon_map():
    assert_equal()

def test_get_uncorrected_no_nt_distances():
    assert_equal(len(get_distances(alns_tpls[0], methods[0], fractions[0], data_types[1])), 2)

def test_get_uncorrected_half_nt_distances():
    assert_equal(len(get_distances(alns_tpls[0], methods[0], fractions[1], data_types[1])), 2)

def test_get_uncorrected_all_aa_distances():
    assert_equal(len(get_distances(alns_tpls[0], methods[1], fractions[2], data_types[0])), 2)

def test_get_jc_all_nt_distances():
    assert_equal(len(get_distances(alns_tpls[0], methods[1], fractions[2], data_types[1])), 2)

def test_get_uncorrected_all_nt_distances():
    assert_equal(len(get_distances(alns_tpls[0], methods[0], fractions[2], data_types[1])), 2)

def test_jc_nt_zero():
    assert_equal(jc_correction((10, 0), data_types[1]), (10, 0))

def test_jc_nan():
    assert_equal(jc_correction((10, 'NaN'), data_types[1]), (10, 'NaN'))

def test_jc_nt_min():
    assert_equal(jc_correction(p_distance(simple_seq, simple_seq), data_types[1]), (10, 0))

def test_jc_aa_min():
    assert_equal(jc_correction(p_distance(simple_aa_seq, simple_aa_seq), data_types[0]), (10, 0))

def test_jc_nt_max():
    assert_equal(jc_correction((10, 10), data_types[1]), (10, (3 / 4 * log(1 - 4 / 3 * -10))))

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