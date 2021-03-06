% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speciesTree.R
\name{TransferDend}
\alias{TransferDend}
\title{Change the format of tree}
\usage{
TransferDend(dend, renameCluster = TRUE, cls.groups)
}
\arguments{
\item{dend}{dendrogram object}

\item{renameCluster}{if rename the leafnode}

\item{cls.groups}{factor contains celltype annotation}
}
\value{
list contains dendrogram obj, cells per leaf node and transfered groups
}
\description{
Change the format of tree
}
