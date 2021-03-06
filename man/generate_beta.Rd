% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{generate_beta}
\alias{generate_beta}
\title{Generates many pairs of optionally correlated SNP effects}
\usage{
generate_beta(M_s, M_p, D, rho)
}
\arguments{
\item{M_s}{Integer. Number of shared direct effects to simulate per phenotype.}

\item{M_p}{Integer or sequence of two integers. Number of private SNPs to
simulate per phenotype pair.}

\item{D}{Integer >= 2. Number of phenotypes to simulate.}

\item{rho}{Float in [-1, 1]. Correlation of each pair of phenotypes.}
}
\description{
Generates many pairs of optionally correlated SNP effects
}
