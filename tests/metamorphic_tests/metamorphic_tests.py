#! /usr/bin/env python3

from nose.tools import assert_true, assert_false, assert_equal, assert_raises

# the data is a simulated nucleotide alignment with 50 taxa, 500 loci each 500 nt long 
# with every taxon having 5 randomly selected loci complemented
# making it easy to flag the outliers and should result in 1% seq removed

# metamorphic tests expect that the results will be the same if
# 1) all sequences are complemented
# 2) order of loci in alignment is scrambled
# 3) order of taxa in alignment is scrambled 