#! /usr/bin/env python3

# coding: utf-8
import random, string
from sys import argv 

from spruceup import aln_parsing, aln_writing


def add_letters(aln_dict):
    for taxon in aln_dict.keys():
        new_taxon = f'{random.choice(string.ascii_uppercase)}{taxon}'
        aln_dict[new_taxon] = aln_dict.pop(taxon)
    return aln_dict


def wrapper(aln_tuple):
    aln_name, aln_dict = aln_tuple
    scram_aln = add_letters(aln_dict)
    aln_writing.write_alignment_file(
        scram_aln, 'fasta', f'scrambled-taxa-{aln_name}', 'dna'
        )


def get_parsed_aln(aln_filename, file_format):
    aln_tuple = aln_parsing.parse_alignment(aln_filename, file_format)
    return aln_tuple


script, aln_filename, file_format = argv
aln_tuple = get_parsed_aln(aln_filename, file_format)
wrapper(aln_tuple)
