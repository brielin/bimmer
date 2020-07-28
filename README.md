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

- R (> 3.5.0)
- [inspre](https://github.com/brielin/inspre)
- dyplyr
- plyr
- purrr
- tidyr
- tibble
- readr
- foreach
- MendelianRandomization
- igraph

To run the tests, simulations and make the corresponding plots, some additional
packages are required:

- mvtnorm
- scales
- egg
- GGally
- huge
- network
- glmnet

Please note that parsing a large number of genome-wide summary statistics files
in R can take substantial time and memory. We recommend at least 64GB and 8
cores if you wish to process 100s of phenotypes.

## Installation

At the moment, the package must be installed with
`devtools::install_github("brielin/bimmer")`. We aim to make the package
available on CRAN eventually.

If you would like to install and run tests,
```
> devtools::install_github("brielin/bimmer", INSTALL_opts="--install-tests")
> install.packages("mvtnorm")  # If not already installed.
> library(testthat)
> library(bimmer)
> test_package("bimmer")
```

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

# Usage Instructions

## Data pre-processing considerations
Bimmer takes GWAS summary statistics as input and currently supports the Neale
lab summary statistics format. The GWAS summary files should
be tab separated (`.tsv`) files and must contain the following columns:

- `minor_AF`: the in-sample or reference minor allele frequency
- `tstat`: the t-statistic associated with the test, sometimes called "Z-score"
- `n_complete_samples`: the total number of individuals for that SNP-trait
  association, total case + control for binary phenotypes
- `beta`: the regression effect size
- `se`: the regression standard error
- `pval`: the p-value associated with the test
- (optional) `low_confidence_variant`: variants marked "TRUE" here
  will be filtered out when considering analyses involving this phenotype
  
Any other columns in the file will be ignored.
  
The main function for importing your data is `bimmer::read_sumstats_neale()`. It
takes as input either 1) a glob file pattern pointing to a set of files on disk
(argument `file_pattern`) or 2) a list of filenames as an R list object
(argument `file_list`). It returns a list of two matrices corresponding to
the harmonized effect size and standard error (`$beta_hat` and `$se_hat`).
*Please note that this function can consume huge memory*. It helps to make sure
in advance that your summary statistics files have the *exact* same variant IDs.
If they do not, this function will raise a warning and may take even longer to
finish.
  
*bimmer requires two sets of summary statistics*. The first set is used to
select SNPs for inclusion in the model fitting, and the second set is used
to fit the model itself. There is *no statistically valid option* for a single
set of summary statistics: it is imperitive that you have two sets.

*bimmer requires LD-pruned summary statistics*. If there are variants in LD in
your summary statistics files, you have two options. 1) remove them using the
method of your choice, or 2) generate a set of "snps_to_use" files. The latter
files (one for each phenotype) should contain new-line separated names of
SNP ids that are not in LD for each phenotype. Then, you can use the function
`bimmer::read_snp_list()` to parse them into R. This function takes as argument
a glob-file pattern and a returns a list of lists of the valid SNPs for each
phenotype.

## Fitting the TCE matrix

The next step is to select instruments using the function
`bimmer::select_snps()`. The primary argument is `sumstats`, simply the output
of `bimmer::read_sumstats_neale()` on the *first* set of summary statistics.
If your sumstats contain variants in LD, you
can use the argument `snps_to_use` and supply it the outpt of
`bimmer::read_snp_list()`. The argument `p_thresh` can be used to set the
p-value cutoff for SNP association. Reasonable values range from 10^-5 to
10^-8, the default is 5*10^-6. The remaining arguments are primarily for
comparison in simulation.

Next you can use the function `bimmer::fit_tce()` to calculate the TCE of every
phenotype on every other. The argument `sumstats` should be supplied with the
output of `bimmer::read_sumstats_neale()` on the *second* set of summary
statistics. The argument `selected_snps` should be supplied with the output of
`bimmer::select_snps()` from the previous step. The argument `mr_method` allows
you to choose the method that you want to use for caculating the TCE. The
default is `egger_w` for weighted Egger regression as described in the
manuscript. The currently implemented options are

- `egger_w` is the recommended option for speed and robustness to
  correlated pleiotropy
- `egger` or `egger_p` for standard or penalized Egger regression if you are
  *not* concerned about correlated pleiotropy
- `mbe` for the "mode best estimator"
- `ps`, `aps`, or `raps` for various versions of the "robust adjusted profile 
  score" which are *very* powerful but *very* susceptible to false positives due
  to correlated pleiotropy.
  
Additional methods for TCE calculation can be added by using a wrapper function
as demonstrated in the body of `bimmer::fit_tce()`. If you have one you like,
feel free to submit a pull request implementing the functionality.

If you are analyzing many phenotypes it is optional but recommended for you to
do some filtering of the TCE matrix using the function `bimmer::filter_tce()`.
This caps the empirical value of the TCE to 1.0 and removes phenotypes that
appear to be problematic as either exposures or outcomes.

## Finding the DCEs

Finally, we convert the TCE matrix into a DCE matrix. If you would like to
use weights (recommended), you must first construct a weight set using
`inspre::make_weights()`. The argument `SE` should be the standard error matrix
(`$SE_tce`) from the output of the function `bimmer::fit_tce()` (or 
`bimmer::filter_tce()` if you chose to use it). The argument `max_min_ratio`
can be used if you have some TCEs with very small SEs to prevent over-fitting
to a handful of data points. The value used in the manuscript of
`max_min_ratio=10000` is a reasonable place to start, but *you may have to
adjust this* if you are having trouble getting a good model fit in the next
step.

The final step is to use `bimmer::fit_inspre()` to fit the DCE matrix. The 
argument `R_tce` should be set to the TCE matrix (`$R_tce`) output by
`bimmer::fit_tce()` (or, again, `bimmer::filter_tce()` if you chose to use it).
The argument `W` can be `NULL`, to use no weights, or set to the output
of `inspre::make_weights()`. If you have a multicore machine (recommended)
you can set `ncores` to speed up the fit but note that this will require a lot
of memory as well.

The arguments `lambda`, `lambda_min_ratio` and `nlambda` control the search for
the regularization strength. You can use `lambda` to either set a single value
or a sequence of values yourself, otherwise it will search over a
logarithmically spaced set of `nlambda` values from 1.0 to `lambda_min_ratio`.
To choose a specific value of lambda for analysis, you can use cross-validation
by supplying the argument `cv_folds` to a non-zero number, we recommend
`cv_folds=10`. If you do this, the result of `bimmer::fit_inspre()` will
contain a vector `$D_hat` with the estimated "stability" for each setting of
`lambda`. We recommend using the value of `lambda` that corresponds to a
`D_hat` of a little bit smaller than 0.025. If all of your estimates for `D_hat`
are above this number, try a sequence with smaller values of `lambda` and vice
versa if all of your estimates are below it. If you have estimates both above
and below 0.025 but none are close to it, try a narrow sequence of lambda values
with less spacing between them.

## Putting it all together

If you put this all together, your R script, `Rmd` or Rstudio command line should
look something like this,

```
sumstats_files_1 <- "/path/to/files/phenotype_*_set_1.tsv"
sumstats_files_2 <- "/path/to/files/phenotype_*_set_2.tsv"
clumped_snps <- "/path/to/files/snps_phenotype_*.txt"

sumstats_1 <- bimmer::read_sumstats_neale(file_pattern = sumstats_files_1)
snps_to_use <- bimmer::read_snp_list(file_pattern = clumped_snps)
selected_snps <- bimmer::select_snps(sumstats = sumstats_1, snps_to_use = snps_to_use)

sumstats_2 <- bimmer::read_sumstats_neale(file_pattern = sumstats_files_2)
tce_result <- bimmer::fit_tce(sumstats = sumstats_2, selected_snps = selected_snps)
tce_filtered <- bimmer::filter_tce(tce_result$R_tce, tce_result$SE_tce)

weights <- inspre::make_weights(SE = tce_filtered$SE_tce, max_min_ratio = 10000)
dce_result <- bimmer::fit_inspre(R_tce = tce_filtered$R_tce, W = weights, cv_folds = 10)
selected_index <- which(cde_res_narrow$D_hat > beta)[1] - 1
lambda <- dce_result$lambda[selected_index]
R_hat <- dce_result$R_hat[ , , selected_index]
U_hat <- dce_result$U[,,selected_index]
```

After this, `R_hat` contains the estimate of the DCEs at the selected index,
`U_hat` contains the shrunken estimates of the TCEs, and `lambda` contains the
regularization parameter chosen by cross-validation.

Because there are many steps in this process, users are highly encouraged to
run each step separately and check the output for reasonable behavior after
every step.



