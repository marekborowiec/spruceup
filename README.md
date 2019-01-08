[![DOI](https://zenodo.org/badge/@@@.svg)](https://zenodo.org/badge/latestdoi/@@@)


```
                    A  
                   TCG  
                  ACGTA  
                    T  
```

# :christmas_tree: seq-spruceup

Tools to discover, visualize, and remove outlier sequences in large multiple sequence alignments. 

If you are using this program, please cite: [this publication](link):
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

# Interface
Once you successfully installed `XX` you will need 1) an alignment in `FASTA`, `PHYLIP` or `NEXUS` format and 2) a configuration file to run it. To run the program from the command line you can type:
```bash
XX my-configuration-file.conf
```
Directory `Examples` contains a template configuration file. It has the following fields:
### [input]
The input category defines parameters of the input alignment and its type.
input_file_name: 
input_format:
data_type:
cores:
method:
distance_method:
window_size:
stride:
fraction:
criterion:
cutoffs:
taxon:
manual_cutoffs:
You can inactivate this option by commenting it out (placing `#` before `manual_cutoffs`).
output_file_aln:
output_format:
report:
log:


[analysis]
[output]



## Examples

