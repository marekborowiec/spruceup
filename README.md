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

This script uses [numpy](link), [scipy](link), [matplotlib](link), [psutil](link) and [tqdm](link).

## Installation and requirements

You can download a zipped GitHub repository, clone it if you have `git` installed on your system, or install using [pip](https://pip.pypa.io/en/latest/installing.html) (recommended) from the [Python Package Index](https://pypi.python.org/pypi/XX/):
```bash
pip install XX
```

`XX` requires you have Python version 3.6 or newer. Dependencies should be installed automatically. If your system does not have Python version 3.6 or newer you will need to [download and install it](http://www.python.org/downloads/). On Linux-like systems (including Ubuntu) you can install it from the command line using

```bash
sudo apt-get install python3.6
```




```bash

```

## Interface
Once you successfully installed `XX` you will need 1) an alignment in `FASTA`, `PHYLIP` or `NEXUS` format and 2) a configuration file to run it. To run the program from the command line you can type:
```bash
spruceup my-configuration-file.conf
```
Directory `Examples` contains a template configuration file. It has the following fields:
### [input]
The `input` category defines parameters of the input alignment and its type.
`input_file_name` is the file path of the alignment to be processed. 
`input_format` indicates which format the alignment file is in. It can be one of five popular formats: `fasta`, `phylip`, `phylip-int` (interleaved PHYLIP), `nexus`, and `nexus-int`. 
`data_type` tells the program whether your alignment contains amino acids (`aa`) or DNA nucleotides (`nt`).  
### [analysis]
The `analysis` category defines parameters used to analyze and clean up your alignment. 
`cores` how many CPU cores to use for distance calculations.
`distance_method` chooses to compute uncorrected distance with `p-distance` or Jukes-Cantor-corrected distance with `jc69`.
`window_size` chooses how many characters (aa/nt) to include in a window in which distances will be calculated.
`stride` chooses how many characters (aa/nt) to move with each slide of the window. Stride of 5 and window size of 20 means that each new window will overlap by 15 characters with the preceding window. Increasing stride will decrease computational burden because
### maybe rename this to overlap and create a function that would do that??!
`fraction` signifies proportion of OTUs/samples that will be used to calculate average distance in each window. With fraction set to `1.0` distances for each OTU will be calculated against all other OTUs in the alignment. With fraction set to `0.5` distances for each OTU will be calculated against a random sample representing 50% of OTUs in the alignment. This setting may help to speed up calculations in alignments with large numbers of taxa.   
`criterion` chooses how outlier distances will be determined. `lognorm` (recommended) means that a [lognormal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution) will be fitted to your distance data for each OTU and cutoffs will be determined by specifying quantile of observations above which sequence will be considered outliers. If you are using `mean` or `median`, simple multiple of those values computed for each OTU will be considered cutoffs for identifying outliers. 
`cutoffs` specifies multiple values considered as cutoffs. For `lognorm` criterion use fractions of `1`, for example `0.9,0.995` etc. If you are using `mean` or `median` as your criterion, use multiples of those values, for example `5` or `15`.  
`manual_cutoffs` is an optional setting that allows manual modifications to cutoffs for individual OTUs. It may prove useful if only one or a few samples have a significant proportion of poorly aligned sequences, skewing their overall cutoff such that they are not being flagged. If you find that this is case, however, you should probably rather be checking your data and pipeline for errors!
### [output]
The `output` category tells the program how and where to save your analysis results.
`output_file_aln` is the name for your trimmed output alignment(s). The actual name saved on your machine will have a prefix signifying cutoff value used.  
`output_format` file format for the trimmed alignment. Choose from `fasta`, `phylip`, `phylip-int`, `nexus`, or `nexus-int`. 
`report` name of files containing information on which sequences were flagged as outliers. The actual name will have a prefix signifying cutoff value used.
`log` is the name of the log with all analysis screen output.






## Examples

