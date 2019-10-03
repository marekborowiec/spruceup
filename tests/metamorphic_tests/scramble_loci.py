#! /usr/bin/env python3

# coding: utf-8
import os
import random
from sys import argv 

from spruceup import aln_writing

# load in alignments
# for each taxon reverse complement 1 seq in all 5 different loci
# write modified alignments 
wd = os.path.dirname(os.path.realpath(__file__)) 
# reverse complement
def split_seq(dna_string, chunk_size):
    return [dna_string[i:i+int(chunk_size)] 
     for i in range(0, len(dna_string), int(chunk_size))]


def rearrange(seq_list, indices):
    new_seq_list = [seq_list[i] for i in indices]
    return new_seq_list


def scramble_alignment(aln_dict, chunk_size):
    # load in aln_dict
    # split each taxon's seq into 500 pieces
    # get a random list of 500 indices
    scrambled_order = random.sample(range(100), k=100)
    for taxon in aln_dict.keys():
        aln_dict[taxon] = split_seq(aln_dict[taxon], chunk_size)
    for taxon in aln_dict.keys():
        aln_dict[taxon] = ''.join(rearrange(aln_dict[taxon], scrambled_order))
    return aln_dict


def scramble_wrapper(aln_tuple, chunk_size):
    aln_name, aln_dict = aln_tuple
    f_name = os.path.basename(aln_name)
    new_aln_name = f'{wd}/scrambled-loci-{f_name}'
    scram_aln = scramble_alignment(aln_dict, chunk_size)
    aln_writing.write_alignment_file(
        scram_aln, 'fasta', new_aln_name, 'dna'
        )
