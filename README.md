[![DOI](https://joss.theoj.org/papers/10.21105/joss.01635/status.svg)](https://doi.org/10.21105/joss.01635)

# spruceup
```
        A  
       TCG  
      ACGTA  
        T 
```

Tools to discover, visualize, and remove outlier sequences in large multiple sequence alignments.

`spruceup` is a Python tool for biologists (bioinformaticians, phylogeneticists, evolutionary biologists) doing inference on phylogenomic sequence alignments. It allows discovery and removal of individual poorly aligned sequences or sequence fragments (alignment rows), which is different from the problem of poorly aligned sequence blocks (alignment columns) commonly addressed by alignment trimming software.

If you are using this program, please cite [this publication](https://joss.theoj.org/papers/10.21105/joss.01635):
```
Borowiec, M.L. (2019) Spruceup: fast and flexible identification, visualization, and removal of outliers from large multiple sequence alignments. Journal of Open Source Software, 4(42), 1635, https://doi.org/10.21105/joss.01635
```

This script uses [numpy](http://www.numpy.org/), [scipy](https://www.scipy.org/index.html), [matplotlib](https://matplotlib.org/), [psutil](https://pypi.org/project/psutil/), [tqdm](https://pypi.org/project/tqdm/), and [treeswift](https://pypi.org/project/treeswift/).

## Installation and requirements

You can download a zipped GitHub repository, clone it if you have `git` installed on your system, or install using [pip](https://pip.pypa.io/en/latest/installing.html) from the [Python Package Index](https://pypi.python.org/pypi/spruceup). Cloning the repository will also download example files used in the tutorial below. 


`spruceup` has been developed and tested under Ubuntu Linux. It requires you have Python version 3.6 or newer. You can use `spruceup` with early releases of Python 3.7, but as of today (12 August 2019), you may experience trouble installing it if using Python 3.7.4. Dependencies should be installed automatically. If your system does not have Python version 3.6 or newer you will need to [download and install it](http://www.python.org/downloads/). On Linux-like systems (including Ubuntu) you can install it from the command line using:

```bash
sudo apt-get install python3.6
```

It is highly recommended that you install `spruceup` in a new Python environment. If you are using the [conda package manager](https://www.anaconda.com/distribution/), type:
```bash
conda create --name spruceup python=3.6
```
And follow instuctions on the prompt.

Now activate your environment (with `conda` you would type `conda activate spruceup`) and install. If you want to download example files, simply download and unzip this repository or, if you have git, clone it with:
```bash
git clone https://github.com/marekborowiec/spruceup.git
```

You can then install `spruceup` from within this cloned repository:
```bash
python setup.py install
```

`spruceup` is also available on [Python Package Index](https://pypi.org/) and you can install it via `pip` without example files:
```bash
pip install spruceup
```

If you do not wish to install `conda`, you can set up a new environment using a different package manager, such as `venv`:
```bash
python3 -m venv spruceup
```

## Interface
Once you successfully installed `spruceup` you will need 1) an alignment in `FASTA`, `PHYLIP` or `NEXUS` format, 2) (optional) a guide tree for your alignment in `NEWICK` format, and 3) configuration file to run the program. To run the program from the command line you can type:
```bash
python -m spruceup my-configuration-file.conf
```
Directory `examples` contains a template configuration file (`/examples/config_example.conf`). This file should be a template for your own analyses. It has the following fields:

### [input]
The `input` category defines parameters of the input alignment and its type.

`input_file_name` is the file path of the alignment to be processed. This can be provided as a file name, if the config is in the same directory as alignment and `spruceup` command will be ran from there as well. Alternatively, you can provide the absolute path. *__This should be a relatively large, concatenated alignment__*. The `spruceup` algorithm works the better the more data it has, so if you have a phylogenomic dataset consisting of many single-locus alignments, you should first concatenate them. You can then split the processed/trimmed alignment back into single-locus alignments using another utility such as [AMAS](https://doi.org/10.7717/peerj.1660).

`input_format` indicates which format the alignment file is in. It can be one of five popular formats: `fasta`, `phylip`, `phylip-int` (interleaved PHYLIP), `nexus`, and `nexus-int`.

`data_type` tells the program whether your alignment contains amino acids (`aa`) or DNA nucleotides (`nt`).

`distances_object_file` is the file name of an existing distance object. By default this saves to the directory where you run the `spruceup` command and so it is better to provide absolute path here. Because sometimes you will want to adjust cutoffs or cutoff criterion and computing distances is the most time-consuming part of the analysis, `spruceup` saves a `json` format file with distances from each analysis. By default this is blank, but if you do have a distance file from a previous analysis and you want to trim your alignment with new cutoffs, supply the `json` file name here. `spruceup` will then run with new trimming cutoffs and/or criterion but without the need to re-calculate distances.

`guide_tree` is a phylogram or cladogram `NEWICK` format file to be used as a guide tree. File name or absolute path. The tree can be inferred using any method and can be fully resolved or contain polytomies. If you do not supply a guide tree the program will still run but without phylogeny it will have less information to identify misaligned sequences. This is particularly important in vartiable algnments with distantly related samples (individual sequences representing individual, taxon, OTU, etc.), where it is more difficult to distinguish genuinely variable sequences from misaligned fragments. In cases where you suspect and are mainly concerned with samples with spuriously long terminal branches it is adivsable to supply topology-only (cladogram) guide tree and/or run the program without a guide tree. 

### [analysis]
The `analysis` category defines parameters used to analyze and clean up your alignment.

`cores` how many CPU cores to use for distance calculations.

`distance_method` chooses to compute uncorrected p-distance with `uncorrected` or Jukes-Cantor-corrected distance with `jc`.

`window_size` chooses how many characters (aa/nt) to include in a window in which distances will be calculated. Default value that works well for most alignments is `20`. 

`overlap` indicates how many characters (aa/nt) each sliding window will be overlapping with preceding window. Overlap of `15` and window size of `20` means that each new window will move 5 positions down the alignment and overlap by 15 characters with the preceding window. Decreasing the overlap will decrease computational burden because fewer windows will be created. Default value is `15` (two thirds of default window size of `20`) but you may want to go lower, to half of window size or even `0` (non-overlapping windows) if your alignment is very large and you want to decrease compute time and memory usage and don't mind sacrificing some precision. 

`fraction` signifies proportion of samples that will be used to calculate average distance in each window. With fraction set to `1.0` distances for each sample will be calculated against all other samples in the alignment. With fraction set to `0.5` distances for each sample will be calculated against a random draw representing 50% of samples in the alignment. Lowering this number will help to speed up calculations in alignments with large numbers of taxa.

`criterion` chooses how outlier distances will be determined. `lognorm` means that a [lognormal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution) will be fitted to your distance data for each sample and cutoffs will be determined by specifying quantile of observations above which sequence will be considered outliers. If you are using `mean`, simple multiple of those values computed for each sample will be considered cutoffs for identifying outliers.

`cutoffs` specifies multiple values considered as cutoffs. If using `lognorm` criterion, use fractions of `1`, for example `0.9,0.995` etc. Default values that should work for most alignments are `0.9,0.95,0.97,0.99,0.995,0.999`. If you are using `mean` as your criterion, use multiples of those values, for example `5,30` etc. Defaults of `5,10,15,18,20,25,30` for `mean` criterion should be useful starting points for most datasets. You can always trim with additional criteria and cutoffs after the initial analysis (see point 5 below under Interpreting the output). If the alignment contains many saturated or poorly aligned sites, a low setting may result in huge amount of data being trimmed from the original alignment. This is time-consuming and you may want to trim your alignment with a more stringent 'block' method before using `spruceup` or remove lower cutoff values from the list.

`manual_cutoffs` is an optional setting that allows manual modifications to cutoffs for individual samples. It may prove useful if only one or a few samples have a significant proportion of poorly aligned sequences, skewing their overall cutoff such that they are not being flagged. If you find that this is case, however, you should probably rather be checking your data and pipeline for errors!

### [output]
The `output` category tells the program how and where to save your analysis results. *Remember that by default any file will be written to directory from which you executed `spruceup` command. To change this behavior provide full paths.* 

`output_file_aln` is the name for your trimmed output alignment(s). The actual name saved on your machine will have a prefix signifying cutoff value used.

`output_format` file format for the trimmed alignment. Choose from `fasta`, `phylip`, `phylip-int`, `nexus`, or `nexus-int`.

`report` name of files containing information on which sequences were flagged as outliers. The actual name will have a prefix signifying cutoff value used.

`log` is the name of the log with all analysis screen output.

## Examples and interpretation of results

To use `spruceup` you will need to run the `spruceup` script from the command line and provide the name of your configuration file as the argument. If you downloaded and installed `spruceup` from source, the example alignment and config files can be found in `examples` subdirectory:
```bash
cd ./examples
python -m spruceup config_example.conf
```
The example alignment is a subset of empirical, anonymized [ultraconseved element or UCE](https://www.ultraconserved.org/) data set generated from insects.   

Once you run the script, the sequence alignment will be divided into a number of windows of the size and overlap you specified. The script will then compute distances for each sample (alignment row: the sequence representing an individual, taxon, OTU etc.) in each window. This is done all-by-all by default or all-by-fraction of samples, if specified. You will see some messages along the way, including a progress bar that will display the number of iterations (windows) and remaining time, as we as the maximum amount of memory used for distance calculation.

![progress-bar](./README_files/progress-bar.png) 
 
Once all distances are calculated, `criterion` and `cutoffs` settings will determine which windows are considered outliers and should be trimmed out of the alignment. When using the `lognorm` criterion, specifying a quantile of `0.99` means that any sequence window that lies above 99th percentile of distances of a given sample to other sample in that window will be deemed an outlier and should be removed. In theory, setting of `0.99` should mean that 1% of all sequence windows will be removed from each sample. In practice, this is not true because real-life sequence data does not perfectly fit into lognormal distribution. Cutoff being constant, certain samples may have many outlier (misaligned) sequence fragments and more than 1% of sequence data removed, while others may have no misaligned fragments and no outliers.

You can now go back to your configuration file and try other cutoffs or methods without the need to re-calculate distances (unless you would like to use different correction or scaling). Simple load the generated `json` file (see below) with the `distances_object_file` option under `[input]` category.

### Interpreting output

`spruceup` produces several types of output:

1. Report files ending with suffix `-report.txt`, one of which is written for each cutoff specified, which indicated by the prefix (e.g `0.95`, `0.99` and so forth). These files contain the distance cutoff value for each sample and which sequence windows were determined to be outliers and removed.

2. Trimmed alignment files ending with suffix `-trimmed.fas/phylip/nexus`, again, one for each cutoff value. These are the alignments with outlier windows removed.

3. Distance distribution `png` plots, one for each sample. These images can be used to examine the distribution of distances for each taxon, its fit to lognormal distribution, and cutoff value placed on each sample given a cutoff.

An example of distances plot is below. The header is the name of the sample. The x-axis indicates distance to other samples, ranging from `0` to the maximum distance that was found for the sample. The y-axis specifies relative number of windows. Blue bars comprise the histogram of distances. Orange line is the fitted lognormal distribution (only shown when using the `lognorm` criterion) and the vertical dashed line indicates the cutoff above which any window will be deemed an outlier and removed.

The first example plot below shows an sample with relatively smooth distance distribution and few sequence windows with extreme values, none of which are greater than `0.25`. You may not be able to see individual window distances as visible histogram bars since the distributions comprise thousands (in this small example) to hundreds of thousands of distances.

![example-plot-good](./README_files/example-plot-good.png)

The following plot shows a sample with less smooth distribution, overall sequences fewer sequences due to missing data, as indicated by the heigth of the histogram bars, and many outlier windows.

![example-plot-poor](./README_files/example-plot-poor.png)

Both examples have been processed under the same cutoffs values, `0.95,0.97,0.99` quantiles of fitted lognormal distribution but in the poorly aligned data the last falls outside of computed distance values. These plots, combined with visual examinaton of report files and alignments should serve you as a guide on what criterion and cutoff values make most sense for your dataset.

4. Log file 

This is a log file that will contain the same information that appears on the terminal screen, excluding progress bar. It will log what  steps were taken by the program, their timing, show information about how many/what proportion of sites were trimmed at each cutoff and output file names that were written. 

5. Distances Python object file

Calculating distances with `spruceup` is often the most time- and memory-consuming part of the process, although trimming very large numbers of positions from the alignment can also take a long time. Because of this `spruceup` writes a `json` format file each time you run an analysis from scratch, allowing you to load it up later and trim with different criterion or cutoff values. Note that distances will be specific for each window size, overlap, and taxon fraction and you will need to re-run the whole analysis if you want to adjust these parameters. Note that the `json` file can be quite large at >150MB per 100,000 windows and 100 taxa.

## Testing

Tests written for `spruceup` code use the Python's standard library module `unittest` and are integrated with `setuptools`. This means that if you downloaded the `spruceup` source code, you can run tests from the top `spruceup` directory after installing with:
```
python setup.py install
python setup.py test
```
These tests are likely to take several minutes because they include metamorphic tests (cf. [Giannoulatou et al. 2014](https://doi.org/10.1186/1471-2105-15-S16-S15)) which involve running `spruceup` on simulated data with known properties. You can speed things up by increasing the number of cores in the metamorphic test configuration file found under `tests/metamorphic_tests/simulation.conf`.
If you want just a quick check of basic functions, run unit tests without the time-consuming metamorphic tests. To do this, from the top `spruceup` directory run:
```
python -m unittest tests.test_simple
```

## Issues and development

If you encounter bugs or problems running the code, [create a new issue on GitHub](https://help.github.com/en/articles/creating-an-issue). Everyone is encouraged to [contribute](https://github.com/firstcontributions/first-contributions).
