[![DOI](https://zenodo.org/badge/@@@.svg)](https://zenodo.org/badge/latestdoi/@@@)

# seq-spruceup
```
        A  
       TCG  
      ACGTA  
        T 
```

Tools to discover, visualize, and remove outlier sequences in large multiple sequence alignments. 

If you are using this program, please cite [this publication](link):
```
```

This script uses [numpy](http://www.numpy.org/), [scipy](https://www.scipy.org/index.html), [matplotlib](https://matplotlib.org/), [psutil](https://pypi.org/project/psutil/) and [tqdm](https://pypi.org/project/tqdm/).

## Installation and requirements

You can download a zipped GitHub repository, clone it if you have `git` installed on your system, or install using [pip](https://pip.pypa.io/en/latest/installing.html) (recommended) from the [Python Package Index](https://pypi.python.org/pypi/seq-spruceup/) (to-do):
```bash
pip install seq-spruceup
```

`seq-spruceup` requires you have Python version 3.6 or newer. Dependencies should be installed automatically. If your system does not have Python version 3.6 or newer you will need to [download and install it](http://www.python.org/downloads/). On Linux-like systems (including Ubuntu) you can install it from the command line using

```bash
sudo apt-get install python3.6
```


```bash

```

## Interface
Once you successfully installed `seq-spruceup` you will need 1) an alignment in `FASTA`, `PHYLIP` or `NEXUS` format and 2) a configuration file to run it. To run the program from the command line you can type:
```bash
spruceup my-configuration-file.conf
```
Directory `examples` contains a template configuration file. It has the following fields:

### [input]
The `input` category defines parameters of the input alignment and its type.

`input_file_name` is the file path of the alignment to be processed.

`input_format` indicates which format the alignment file is in. It can be one of five popular formats: `fasta`, `phylip`, `phylip-int` (interleaved PHYLIP), `nexus`, and `nexus-int`.

`data_type` tells the program whether your alignment contains amino acids (`aa`) or DNA nucleotides (`nt`).

`distances_object_file` is the file name of an existing distance object. Because sometimes you will want to adjust cutoffs or cutoff criterion and computing distances is the most time-consuming part of the analysis, `spruceup` saves a `json` format file with distances from each analysis. By default this is blank, but if you do have a distance file from a previous analysis and you want to trim your alignment with new cutoffs, supply the `json` file name here. `spruceup` will then run with new trimming cutoffs and/or criterion but without the need to re-calculate distances.

### [analysis]
The `analysis` category defines parameters used to analyze and clean up your alignment.

`cores` how many CPU cores to use for distance calculations.

`distance_method` chooses to compute uncorrected distance with `p-distance` or Jukes-Cantor-corrected distance with `jc69`.

`window_size` chooses how many characters (aa/nt) to include in a window in which distances will be calculated. Default value that works well for most alignments is `20`. 

`overlap` indicates how many characters (aa/nt) each sliding window will be overlapping with preceding window. Stride of 15 and window size of 20 means that each new window will move 5 positions down the alignment and overlap by 15 characters with the preceding window. Increasing overlap will decrease computational burden because fewer windows will be created. Default value is `15` (two thirds of default window size of `20`) but you may want to go lower, to half of window size or even `0` (non-overlapping windows) if your alignment is very large and you want to decrease compute time and memory usage and don't mind sacrificing some precision. 

`fraction` signifies proportion of OTUs/samples that will be used to calculate average distance in each window. With fraction set to `1.0` distances for each OTU will be calculated against all other OTUs in the alignment. With fraction set to `0.5` distances for each OTU will be calculated against a random sample representing 50% of OTUs in the alignment. Lowering this number will help to speed up calculations in alignments with large numbers of taxa.

`criterion` chooses how outlier distances will be determined. If you are using `mean` (recommended), simple multiple of those values computed for each OTU will be considered cutoffs for identifying outliers. `lognorm` means that a [lognormal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution) will be fitted to your distance data for each OTU and cutoffs will be determined by specifying quantile of observations above which sequence will be considered outliers.

`cutoffs` specifies multiple values considered as cutoffs. If you are using `mean` as your criterion, use multiples of those values, for example `5,30` etc. For `lognorm` criterion use fractions of `1`, for example `0.9,0.995` etc.  Default values that should work for most alignments are `0.9,0.95,0.97,0.99,0.995,0.999` for `lognorm` and `5,10,15,20,30,50` for `mean`. If the alignment contains many saturated or poorly aligned sites, a low setting may result in huge amount of data being trimmed from the original alignment. This is time-consuming and you may want to trim your alignment with a more stringent 'block' method before using `spruceup` or remove lower cutoff values from the list.

`manual_cutoffs` is an optional setting that allows manual modifications to cutoffs for individual OTUs. It may prove useful if only one or a few samples have a significant proportion of poorly aligned sequences, skewing their overall cutoff such that they are not being flagged. If you find that this is case, however, you should probably rather be checking your data and pipeline for errors!

### [output]
The `output` category tells the program how and where to save your analysis results.

`output_file_aln` is the name for your trimmed output alignment(s). The actual name saved on your machine will have a prefix signifying cutoff value used.

`output_format` file format for the trimmed alignment. Choose from `fasta`, `phylip`, `phylip-int`, `nexus`, or `nexus-int`.

`report` name of files containing information on which sequences were flagged as outliers. The actual name will have a prefix signifying cutoff value used.

`log` is the name of the log with all analysis screen output.

## Examples and interpretation of results

To use `seq-spruceup` you will need to run the `spruceup` script from the command line and provide the name of your configuration file as the argument:
```bash
spruceup.py my-configuration-file.conf
```
Once you run the script, your sequence alignment will be divided into a number of windows of the size and overlap you specified. The script will then compute distances for each OTU (sample/taxon) in each window. This is done all-by-all by default or all-by-fraction of OTUs, if specified. You will see some messages along the way, including a progress bar that will display the number of iterations (windows) and remaining time, as we as the maximum amount of memory used for distance calculation.

![progress-bar](progress-bar.png) 
 
Once all distances are calculated, `criterion` and `cutoffs` settings will determine which windows are considered outliers and should be trimmed out of the alignment. When using the `lognorm` criterion, specifying a quantile of `0.99` means that any sequence window that lies above 99th percentile of distances of a given OTU to other OUTs in that window will be deemed an outlier and should be removed. In theory, setting of `0.99` should mean that 1% of all sequence windows will be removed from each OTU. In practice, this is not true because real-life sequence data does not perfectly fit into lognormal distribution. Cutoff being constant, certain OTUs may have many outlier (misaligned) sequence fragments and more than 1% of sequence data removed, while others may have no misaligned fragments and no outliers.

### Interpreting the output

`spruceup` produces several types of output:

1. Report files ending with suffix `-report.txt`, one of which is written for each cutoff specified, which indicated by the prefix (e.g `0.95`, `0.99` and so forth). These files contain the distance cutoff value for each OTU and which sequence windows were determined to be outliers and removed.

2. Trimmed alignment files ending with suffix `-trimmed.fas/phylip/nexus`, again, one for each cutoff value. These are the alignments with outlier windows removed.

3. Distance distribution `png` plots, one for each OTU and cutoff value. These images can be used to examine the distribution of distances for each taxon, its fit to lognormal distribution, and cutoff value placed on each OTU given a cutoff.

An example of distances plot is below. The header is the name of the OTU. The x-axis indicates distance to other OTUs, ranging from `0` to `1`. The x-axis is limited to the maximum distance that was found for the OTU. The y-axis specifies number of windows. Blue bars comprise the histogram of distances. Orange line is the fitted lognormal distribution (only shown when using the `lognorm` criterion) and the vertical dashed line indicates the cutoff above which any window will be deemed an outlier and removed.

The first example plot below shows an OTU with relatively smooth distance distribution and few sequence windows with extreme values, none of which are greater than `0.2`. You may not be able to see individual window distances as visible histogram bars since the distributions comprise thousands (in this small example) to hundreds of thousands of distances.

![example-plot-good](example-plot-good.png)

The following plot shows an OTU with less smooth distribution, overall sequences fewer sequences due to missing data, as indicated by the heigth of the histogram bars, and many outlier windows.

![example-plot-poor](example-plot-poor.png)

Both examples have the same cutoff value, `0.97` quantile of fitted lognormal distribution. These plots, combined with visual examinaton of report files and alignments should serve you as a guide on what criterion and cutoff values make most sense for your dataset.

4. Log file (to-do)

5. Distances Python object file

Calculating distances with `spruceup` is often the most time- and memory-consuming part of the process, although trimming very large numbers of positions from the alignment can also take a long time. Because of this `spruceup` writes a `json` format file each time you run an analysis from scratch, allowing you to load it up later and trim with different criterion or cutoff values. Note that distances will be specific for each window size, overlap, and taxon fraction and you will need to re-run the whole analysis if you want to adjust these parameters.