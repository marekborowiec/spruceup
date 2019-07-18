#! /usr/bin/env python3

# coding: utf-8
import random
from sys import argv 

import numpy as np

from spruceup import aln_parsing, aln_writing

# load in alignments
# for each taxon reverse complement 1 seq in all 5 different loci
# write modified alignments 

def complement(dna_string):
    replacement1 = dna_string.replace('A', 't')
    replacement2 = replacement1.replace('T', 'a')
    replacement3 = replacement2.replace('C', 'g')
    replacement4 = replacement3.replace('G', 'c')
    return replacement4.upper()

# reverse complement
def split_seq(dna_string, locus_size):
    return [dna_string[i:i+int(locus_size)] 
     for i in range(0, len(dna_string), int(locus_size))]


def complement_alignment(aln_dict, locus_size, no_loci_to_complement):
    for taxon in aln_dict.keys():
        aln_dict[taxon] = split_seq(aln_dict[taxon], locus_size)
    for taxon in aln_dict.keys():
        all_loci = aln_dict[taxon]
        loci_to_compl = np.random.choice(range(len(all_loci)), 
                                          int(no_loci_to_complement), 
                                          replace=False)
        for i in loci_to_compl:
            all_loci[i] = complement(all_loci[i])
        aln_dict[taxon] = ''.join(aln_dict[taxon])
    return aln_dict


def complement_wrapper(aln_tuple, locus_size, no_loci_to_complement):
    aln_name, aln_dict = aln_tuple
    scram_aln = complement_alignment(aln_dict, locus_size, no_loci_to_complement)
    aln_writing.write_alignment_file(
        scram_aln, 'fasta', f'complemented-loci-{aln_name}', 'dna'
        )


def get_parsed_aln(aln_filename, file_format):
    aln_tuple = aln_parsing.parse_alignment(aln_filename, file_format)
    return aln_tuple


script, aln_filename, file_format, locus_size, no_loci_to_complement = argv
aln_tuple = get_parsed_aln(aln_filename, file_format)
complement_wrapper(aln_tuple, locus_size, no_loci_to_complement)
