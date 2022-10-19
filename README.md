# Bidirectional mediated Mendelian randomization (bimmer)

Bideirectional mediated Mendelian randomization is a method for estimating
causal networks on complex traits from genome-wide association study
summary statistics. For more information on the method please consult the paper
[Phenome-scale causal network discovery with bidirectional mediated Mendelian randomization](https://www.biorxiv.org/content/10.1101/2020.06.18.160176v2).
This user documentation is a work in progress. Please check back in the future
for a complete guide on how to run this software on your data from
pre-processing through model selection. 


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

The first step in applying *bimmer* is to estimate the phenotype-phenotype
causal effects. You need two matrices, one containing the effect estimates
(we refer to this as `R_tce`) and the other containing standard errors
(we refer to this as `SE_tce`). These are generally asymmetric matrices with
rows indicating the cause (sometimes called exposure) and columns indicating
the effect (sometimes called outcome). The effect estimates need to be
bidirectional, that is, each phenotype should appear as both a cause and
an effect.

In this context, we use genetic instrumental variables
to estimate the effect using bidirectional Mendelian randomization. In theory
you can use any MR method that you want, but it is important to consider the
possibility of correlated horizontal pleiotropy. In the original preprint,
we introduce a method called *Welch-weighted Egger regression* that can 
reduce the effects of correlated pleiotropy in large-scale analyses of
biobank-style data, and that is the method that we recommend using prior
to network estmation. The WWER github is
[https://github.com/brielin/WWER](https://github.com/brielin/WWER). The WWER
package includes several tools for parsing and preprocessing biobank
summary statistics and includes options for using alternative MR methods.
Thus, even if you would prefer to use a different MR method, I recommend
looking at the WWER package for preprocessing your data. For more
information on WWER specifically, see our 
[AJHG publication](https://www.sciencedirect.com/science/article/pii/S0002929721003839).

If you are using published MR results, you will need to write custom code to
parse the data into the pair of matrices described above. Keep in mind the scaling
of the data. *bimmer* assumes that the causal effect estimates are variance-normalized,
that is, the causal effect estimate `R` represents a per-variance effect of
the exposure on the outcome. This may require normalizing the estimate by the
variances of the exposure and outcome.

## Fitting the network

Once you have your normalized `R_tce` and `SE_tce` matrices, you are ready to fit
the network. In some cases, you may consider filtering these matrices for more
stable network estimates. For example, you may filter pathological `R` values
(those much greater or less than +-1) or rows/columns of the data matrix with
many missing entries. We provide a helper function `bimmer::filter_tce()` to
assisst with this, see for example its use in `Rmd/ukbb_analysis.Rmd`.

If you would like to use weights (recommended), you must first construct a weight set using
`inspre::make_weights()`. The argument `SE` should be the standard error matrix
`SE_tce`. The argument `max_min_ratio`
can be used if you have some TCEs with very small SEs to prevent over-fitting
to a handful of data points. The value used in the manuscript of
`max_min_ratio=10000` could be a reasonable place to start, but you may have to
adjust this if you are having trouble getting a good model fit in the next
step.

The final step is to use `bimmer::fit_inspre()` to fit the network.
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
# Given phenotype by phenotype matrixes R_tce, SE_tce

beta = 0.025
weights <- inspre::make_weights(SE = SE_tce, max_min_ratio = 10000)
network_result <- bimmer::fit_inspre(R_tce = R_tce, W = weights, cv_folds = 10)
selected_index <- which(network_result$D_hat > beta)[1] - 1
lambda <- dce_result$lambda[selected_index]
R_hat <- dce_result$R_hat[ , , selected_index]
U_hat <- dce_result$U[,,selected_index]
```

After this, `R_hat` contains the estimate of the network at the selected index,
`U_hat` contains the shrunken estimates of the TCEs, and `lambda` contains the
regularization parameter chosen by cross-validation.

Because there are many steps in this process, users are highly encouraged to
run each step separately and check the output for reasonable behavior after
every step.
