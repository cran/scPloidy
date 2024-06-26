% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATAC-seq.R
\name{ploidy}
\alias{ploidy}
\title{Infer Ploidy from ATAC-seq Fragment Overlap}
\usage{
ploidy(
  fragmentoverlap,
  levels,
  s = 100,
  epsilon = 1e-08,
  subsamplesize = NULL,
  dobayes = FALSE,
  prop = 0.9
)
}
\arguments{
\item{fragmentoverlap}{Frequency of fragment overlap in each cell
computed by the function \code{fragmentoverlapcount}.}

\item{levels}{Possible values of ploidy. For example,
\code{c(2, 4)} if the cells can be diploids or tetraploids.
The values must be larger than one.}

\item{s}{Seed for random numbers used in EM algorithm.}

\item{epsilon}{Convergence criterion for the EM algorithm.}

\item{subsamplesize}{EM algorithm becomes difficult to converge
when the number of cells is very large.
By setting the parameter (e.g. to 1e4),
we can run EM algorithm iteratively,
first for \code{subsamplesize} randomly sampled cells,
next for twice the number of cells in repetition.
The inferred lambda/theta parameters are used as the initial value
in the next repetition.}

\item{dobayes}{(experimental feature) Whether to perform Bayesian inference,
which takes long computation time.}

\item{prop}{Proportion of peaks that can be fitted with binomal
distribution in ploidy.bayes. The rest of peaks are allowed to
have depth larger than the ploidy.}
}
\value{
A data.frame with each row corresponding to a cell.
For each cell, its barcode, ploidy inferred by 1) moment method,
2) the same with additional K-means clustering,
3) EM algorithm of mixture, and, optionally,
4) Bayesian inference are given.
I recommend using \code{ploidy.moment} or \code{ploidy.em}.
When \code{fragmentoverlapcount} was computed with \code{dobptonext=TRUE},
we only use the chromosomal sites with chromatin accessibility in top 10%.
This requires longer computation time.
}
\description{
Infer Ploidy from ATAC-seq Fragment Overlap
}
