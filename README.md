# Bidirectional mediated Mendelian randomization (bimmer)

Bideirectional mediated Mendelian randomization is a method for estimating
causal networks on complex traits from genome-wide association study
summary statistics. For more information on the method please consult the paper
[Phenome-scale causal network discovery with bidirectional mediated Mendelian randomization](https://www.biorxiv.org/content/10.1101/2020.06.18.160176v2).
This user documentation is a work in progress. Please check back in the future
for a complete guide on how to run this software on your data from
pre-processing through model selection. 

## Requirements
bimmer is written in the R programming language and depends on the following
packages:

- R (> 3.6)
- [inspre](https://github.com/brielin/inspre)
- dyplyr
- purrr
- tidyr
- tibble
- readr
- foreach
- MendelianRandomization
- igraph

To run the simulations and make the corresponding plots, some additional
packages are required:

- mvtnorm
- scales
- egg
- GGally
- huge
- network

To run the tests you will also need the testthat package.

## Overview of the repository
- `R` contains `*.R` files with functions for preprocessing and loading data. The
main files are `sumstats.R`, which contains functions for loading and preprocessing
GWAS summary statistics, and `fit-model.R`, which contains functions for fitting
the model to data
- `Rmd` contains the `*.Rmd` files with sections corresponding to the processing
of the UKBB data and simulations reported in the manuscript
- `ukbb_analysis` contains helper scripts used in the analysis of the UK
biobank data for the manuscript
- `tests` contains tests for the functions in `fit-model.R` and `simulate.R`
- `man` contains the function-level documentation for the package
- `manuscript` contains the LaTeX files and figures for the manuscript


