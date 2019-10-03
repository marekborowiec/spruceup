#! /usr/bin/env python3

# coding: utf-8
import os
from sys import argv 

from spruceup import aln_writing
# load in alignments
# for each taxon reverse complement 1 seq in all 5 different loci
# write modified alignments 

wd = os.path.dirname(os.path.realpath(__file__)) 
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
    f_name = os.path.basename(aln_name)
    new_aln_name = f'{wd}/complemented-{f_name}'
    compl_aln = complement_alignment(aln_dict)
    aln_writing.write_alignment_file(
        compl_aln, 'fasta', new_aln_name, 'dna'
        )
