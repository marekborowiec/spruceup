#! /usr/bin/env python3
import glob, os, unittest, subprocess

from spruceup import aln_parsing
from tests.metamorphic_tests import complement, scramble_loci, scramble_taxa
# the data is a simulated nucleotide alignment with 41 taxa, 100 loci each 500 nt long
# with every taxon having 5 randomly selected loci complemented
# making it easy to flag the outliers
# and should result in 20,500 nt or 1% total seq removed

# metamorphic tests expect that the results will be the same if
# 1) all sequences are complemented
# 2) order of loci in alignment is scrambled
# 3) order of taxa in alignment is scrambled

# these tests are likely to take several minutes

wd = os.path.dirname(os.path.realpath(__file__)) 

class MetamorphicTests(unittest.TestCase):
    def setUp(self):
        aln_filename = 'simulated.fas'
        file_format = 'fasta'
        aln_tuple = aln_parsing.parse_alignment(f'{wd}/{aln_filename}', file_format)
        complement.complement_wrapper(aln_tuple) 
        subprocess.run(['cp', f'{wd}/simulation.conf', f'{wd}/complemented-simulation.conf'], cwd=wd)
        subprocess.run(
            [
                'sed',
                '-i',
                's/simulated/complemented-simulated/g',
                f'{wd}/complemented-simulation.conf',
            ],
            cwd=wd
        )

        scrl_aln_filename = f'{wd}/complemented-simulated.fas'
        chunk_size = '500'
        scrl_aln_tuple = aln_parsing.parse_alignment(scrl_aln_filename, file_format)
        scramble_loci.scramble_wrapper(scrl_aln_tuple, chunk_size)
        subprocess.run(
            ['cp', f'{wd}/simulation.conf', f'{wd}/scrambled-loci-complemented-simulation.conf'],
            cwd=wd
        )
        subprocess.run(
            [
                'sed',
                '-i',
                's/simulated/scrambled-loci-complemented-simulated/g',
                'scrambled-loci-complemented-simulation.conf',
            ], 
            cwd=wd
        )

        scrt_aln_filename = f'{wd}/scrambled-loci-complemented-simulated.fas'
        scrt_aln_tuple = aln_parsing.parse_alignment(scrt_aln_filename, file_format)
        scramble_taxa.wrapper(scrt_aln_tuple)
        subprocess.run(
            [
                'cp',
                f'{wd}/simulation.conf',
                f'{wd}/scrambled-taxa-scrambled-loci-complemented-simulation.conf',
            ],
            cwd=wd
        )
        subprocess.run(
            [
                'sed',
                '-i',
                's/simulated/scrambled-taxa-scrambled-loci-complemented-simulated/g',
                f'{wd}/scrambled-taxa-scrambled-loci-complemented-simulation.conf',
            ],
            cwd=wd
        )
        

    def test_correct_removal_baseline(self):
        baseline = subprocess.run(
            ['python', '-m', 'spruceup', 'simulation.conf'],
            cwd=wd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertTrue('Sites removed: 20500' in str(baseline.stderr))
        self.assertTrue('Removed 1.00%' in str(baseline.stderr))


    def test_correct_removal_complemented(self):
        complemented = subprocess.run(
            ['python', '-m', 'spruceup', 'complemented-simulation.conf'],
            cwd=wd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertTrue('Sites removed: 20500' in str(complemented.stderr))
        self.assertTrue('Removed 1.00%' in str(complemented.stderr))


    def test_correct_removal_loci_randomized(self):
        scrambled_loci = subprocess.run(
            ['python', '-m', 'spruceup', 'scrambled-loci-complemented-simulation.conf'],
            cwd=wd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertTrue('Sites removed: 20500' in str(scrambled_loci.stderr))
        self.assertTrue('Removed 1.00%' in str(scrambled_loci.stderr))


    def test_correct_removal_taxa_randomized(self):
        scrambled_taxa = subprocess.run(
            [
                'python', '-m', 'spruceup', 
                'scrambled-taxa-scrambled-loci-complemented-simulation.conf',
            ],
            cwd=wd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertTrue('Sites removed: 20500' in str(scrambled_taxa.stderr))
        self.assertTrue('Removed 1.00%' in str(scrambled_taxa.stderr))


    def remove_wrapper(self, file_list):
        for fn in file_list:
            try:
                os.remove(fn)
            except:
                print(f'Error while deleting file: {fn}')



    def tearDown(self):
        print('Cleaning up directory ...')
        pngs = glob.glob('{}/*.png'.format(wd))
        cnfs = glob.glob('{}/*-simulation.conf'.format(wd))
        outs = glob.glob('{}/*lognorms*'.format(wd))
        logs = glob.glob('{}/*.log'.format(wd))
        jsns = glob.glob('{}/*.json'.format(wd))
        mtcs = glob.glob('{}/*-simulated.fas'.format(wd))
        self.remove_wrapper(pngs)
        self.remove_wrapper(cnfs)
        self.remove_wrapper(outs)
        self.remove_wrapper(logs)
        self.remove_wrapper(jsns)
        self.remove_wrapper(mtcs)