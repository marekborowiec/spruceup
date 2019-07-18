#! /usr/bin/env python3
import glob, os, subprocess

from nose.tools import assert_true

# the data is a simulated nucleotide alignment with 41 taxa, 100 loci each 500 nt long
# with every taxon having 5 randomly selected loci complemented
# making it easy to flag the outliers
# and should result in 20,500 nt or 1% total seq removed

# metamorphic tests expect that the results will be the same if
# 1) all sequences are complemented
# 2) order of loci in alignment is scrambled
# 3) order of taxa in alignment is scrambled

# these tests are likely to take several minutes

def setup_module():
    subprocess.run(['./complement.py', 'simulated.fas', 'fasta'])
    subprocess.run(['cp', 'simulation.conf', 'complemented-simulation.conf'])
    subprocess.run(
        [
            'sed',
            '-i',
            's/simulated/complemented-simulated/g',
            'complemented-simulation.conf',
        ]
    )
    subprocess.run(
        ['./scramble-loci.py', 'complemented-simulated.fas', 'fasta', '500']
    )
    subprocess.run(
        ['cp', 'simulation.conf', 'scrambled-loci-complemented-simulation.conf']
    )
    subprocess.run(
        [
            'sed',
            '-i',
            's/simulated/scrambled-loci-complemented-simulated/g',
            'scrambled-loci-complemented-simulation.conf',
        ]
    )
    subprocess.run(
        [
            './scramble-taxa.py',
            'scrambled-loci-complemented-simulated.fas',
            'fasta',
        ]
    )
    subprocess.run(
        [
            'cp',
            'simulation.conf',
            'scrambled-taxa-scrambled-loci-complemented-simulation.conf',
        ]
    )
    subprocess.run(
        [
            'sed',
            '-i',
            's/simulated/scrambled-taxa-scrambled-loci-complemented-simulated/g',
            'scrambled-taxa-scrambled-loci-complemented-simulation.conf',
        ]
    )


def teardown_module():
    print('Cleaning up directory ...')
    pngs = glob.glob('*.png')
    cnfs = glob.glob('*-simulation.conf')
    outs = glob.glob('*lognorms*')
    logs = glob.glob('*.log')
    jsns = glob.glob('*.json')
    mtcs = glob.glob('*-simulated.fas')
    remove_wrapper(pngs)
    remove_wrapper(cnfs)
    remove_wrapper(outs)
    remove_wrapper(logs)
    remove_wrapper(jsns)
    remove_wrapper(mtcs)
    

def test_correct_removal_baseline():
    baseline = subprocess.run(
        ['spruceup.py', 'simulation.conf'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert_true('Sites removed: 20500' in str(baseline.stderr))
    assert_true('Removed 1.00%' in str(baseline.stderr))


def test_correct_removal_complemented():
    complemented = subprocess.run(
        ['spruceup.py', 'complemented-simulation.conf'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert_true('Sites removed: 20500' in str(complemented.stderr))
    assert_true('Removed 1.00%' in str(complemented.stderr))


def test_correct_removal_loci_randomized():
    scrambled_loci = subprocess.run(
        ['spruceup.py', 'scrambled-loci-complemented-simulation.conf'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert_true('Sites removed: 20500' in str(scrambled_loci.stderr))
    assert_true('Removed 1.00%' in str(scrambled_loci.stderr))


def test_correct_removal_taxa_randomized():
    scrambled_taxa = subprocess.run(
        [
            'spruceup.py',
            'scrambled-taxa-scrambled-loci-complemented-simulation.conf',
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert_true('Sites removed: 20500' in str(scrambled_taxa.stderr))
    assert_true('Removed 1.00%' in str(scrambled_taxa.stderr))


def remove_wrapper(file_list):
    for fn in file_list:
        try:
            os.remove(fn)
        except:
            print(f'Error while deleting file: {fn}')
