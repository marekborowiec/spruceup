#! /usr/bin/env python
# -*- coding: utf-8 -*-
import random

import dendropy
from dendropy.interop import seqgen
import numpy as np

POSSIBLE_STATE_FREQS = [0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
POSSIBLE_GENERAL_RATES = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
POSSIBLE_SCALE_BRANCH_LENS = (0.1, 50)
SEQ_LENGTH = 500

def simulate_matrix_wrapper(num):
    tree = dendropy.Tree.get(
        path="for-simulation.newick",
        schema="newick")
    for i in range(num):
        freqs = sum_to_one(random.choices(POSSIBLE_STATE_FREQS, k=4))
        rates = random.choices(POSSIBLE_GENERAL_RATES, k=6)
        scale = np.random.exponential(scale=1)
        #scale = random.uniform(*POSSIBLE_SCALE_BRANCH_LENS)
        seq_length = SEQ_LENGTH
        print(freqs, rates, scale)
        fasta_string = simulate_gtr_matrix(tree, seq_length, freqs, rates, scale)
        fn = f'sim{i}.fas'
        with open(fn, 'w') as f:
            f.write(fasta_string)


def simulate_gtr_matrix(tree, seq_length, frequencies, rates, branch_scale):
    s = seqgen.SeqGen()
    s.char_model = seqgen.SeqGen.GTR
    s.state_freqs = frequencies
    s.general_rates = rates
    s.scale_branch_lens = branch_scale
    s.seq_len = seq_length
    d = s.generate(tree)
    fasta_string = d.char_matrices[0].as_string('fasta')
    return fasta_string


def sum_to_one(random_list):
    summing_to_one = [round((e / sum(random_list)), 2) for e in random_list]
    return summing_to_one


simulate_matrix_wrapper(100)

