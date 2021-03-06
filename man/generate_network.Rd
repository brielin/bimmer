% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{generate_network}
\alias{generate_network}
\title{Generates and optionally normalizes causal graph.}
\usage{
generate_network(
  D,
  graph = "random",
  orient = "random",
  prob = 2/D,
  g = D/20,
  v = 0.2
)
}
\arguments{
\item{D}{Integer. Number of phenotypes to simulate.}

\item{graph}{String, one of "random", "hub" or "scale-free". Specifies the
kind of graph to be generated.}

\item{orient}{String, one of "random", "towards" or "away". Specifies the
edge orientation strategy. Randomly, towards high degree nodes, or away
from high degree nodes.}

\item{prob}{Float between 0 and 1. Specifies the probability of edge
inclusion for "random" graphs. Default is 2/D.}

\item{g}{Integer between 1 and D. Specifies the number of hubs for hub
graphs. Default is D/20.}

\item{v}{Float between 0 and 1. Roughly corresponds to edge weight, see
huge.generator for more information. Default 0.3.}
}
\description{
This is a wrapper around the huge.generator function from that package. After
generating the adjacency matrix, this orients the sign randomly and then the
edges according to the specified method.
}
