#! /usr/bin/env python3

# coding: utf-8
from sys import argv 

from spruceup import aln_parsing, aln_writing

# load in alignments
# for each taxon reverse complement 1 seq in all 5 different loci
# write modified alignments 


# reverse complement
def complement(dna_string):
    replacement1 = dna_string.replace('A', 't')
    replacement2 = replacement1.replace('T', 'a')
    replacement3 = replacement2.replace('C', 'g')
    replacement4 = replacement3.replace('G', 'c')
    return replacement4.upper()


def complement_alignment(aln_dict):
    for taxon in aln_dict.keys():
        aln_dict[taxon] = complement(aln_dict[taxon])
    return aln_dict


def complement_wrapper(aln_tuple):
    aln_name, aln_dict = aln_tuple
    compl_aln = complement_alignment(aln_dict)
    aln_writing.write_alignment_file(
        compl_aln, 'fasta', f'complemented-{aln_name}', 'dna'
        )


def get_parsed_aln(aln_filename, file_format):
    aln_tuple = aln_parsing.parse_alignment(aln_filename, file_format)
    return aln_tuple


script, aln_filename, file_format = argv
aln_tuple = get_parsed_aln(aln_filename, file_format)
complement_wrapper(aln_tuple)
