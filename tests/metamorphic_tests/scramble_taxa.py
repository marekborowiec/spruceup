#! /usr/bin/env python3

# coding: utf-8
import os
import random, string
from sys import argv 

from spruceup import aln_writing

wd = os.path.dirname(os.path.realpath(__file__)) 
def add_letters(aln_dict):
    for taxon in aln_dict.keys():
        new_taxon = f'{random.choice(string.ascii_uppercase)}{taxon}'
        aln_dict[new_taxon] = aln_dict.pop(taxon)
    return aln_dict


def wrapper(aln_tuple):
    aln_name, aln_dict = aln_tuple
    f_name = os.path.basename(aln_name)
    new_aln_name = f'{wd}/scrambled-taxa-{f_name}'
    scram_aln = add_letters(aln_dict)
    aln_writing.write_alignment_file(
        scram_aln, 'fasta', new_aln_name, 'dna'
        )
