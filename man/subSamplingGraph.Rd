% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speciesTree.R
\name{subSamplingGraph}
\alias{subSamplingGraph}
\title{subsampling graph}
\usage{
subSamplingGraph(
  g,
  method = rleiden.detection,
  stability.subsamples = 10,
  stability.subsampling.fraction = 0.95,
  saveGraph = T,
  prefix = NULL
)
}
\arguments{
\item{g}{igraph obj}

\item{method}{clustering algorithm (e.g. walktrap, leiden etc)}

\item{saveGraph}{if save the susample result}

\item{stability.subsample}{# of subsampling}
}
\description{
subsampling graph
}
