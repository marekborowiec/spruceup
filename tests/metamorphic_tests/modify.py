#! /usr/bin/env python3

# coding: utf-8
import glob, random

import aln_parsing, aln_writing

# load in alignments
# for each taxon reverse complement 1 seq in 5 different loci
# write modified alignments 
ALN_FILES =  glob.glob('sim*.fas')
FORMAT = 'fasta'
TAXA = ['T1', 'T10', 'T11', 'T2', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 
'T27', 'T28', 'T29', 'T3', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 
'T37', 'T38', 'T39', 'T4', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 
'T47', 'T48', 'T49', 'T5', 'T50', 'T6', 'T7', 'T8', 'T9']
NO_SEQ_TO_MODIFY = 500

print(len(TAXA))
print(len(ALN_FILES))

# reverse complement
def complement(dna_string):
    replacement1 = dna_string.replace('A', 't')
    replacement2 = replacement1.replace('T', 'a')
    replacement3 = replacement2.replace('C', 'g')
    replacement4 = replacement3.replace('G', 'c')
    return replacement4.upper()


def complement_alignment(aln_dict, taxon):
    aln_dict[taxon] = complement(aln_dict[taxon])
    return aln_dict


def complement_wrapper(aln_tuples, taxa, no_seq_to_modify):
    for taxon in taxa:
        #alns_to_modify = random.choices(aln_tuples, k=no_seq_to_modify, replace=False)
        for aln in aln_tuples:
            aln_name, aln_dict = aln
            compl_aln = complement_alignment(aln_dict, taxon)
            print(f'Scrambling alignment: {aln_name} for taxon {taxon}')
            aln_writing.write_alignment_file(
                compl_aln, 'fasta', f'complemented-{aln_name}', 'dna'
            )
            #aln_tuples.remove(aln)


def get_parsed_alns(aln_filenames, file_format):
    aln_tuples = [aln_parsing.parse_alignment(aln, file_format) for aln in aln_filenames]
    return aln_tuples


def modify(aln_filenames, file_format, taxa, no_seq_to_modify):
    aln_tuples = get_parsed_alns(aln_filenames, file_format)
    complement_wrapper(aln_tuples, taxa, no_seq_to_modify)

modify(ALN_FILES, FORMAT, TAXA, NO_SEQ_TO_MODIFY)

