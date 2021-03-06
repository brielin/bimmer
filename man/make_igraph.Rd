% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit-model.R
\name{make_igraph}
\alias{make_igraph}
\title{Creates an igraph from CDE matrix.}
\usage{
make_igraph(R_cde, min_edge_value = 0.01, max_edge_value = 0.999)
}
\arguments{
\item{R_cde}{Matrix of causal effects.}

\item{min_edge_value}{Minimum edge strength for pruning.}

\item{max_edge_value}{Set edges above this number to this.}
}
\description{
Creates an igraph from CDE matrix.
}
