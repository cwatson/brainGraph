\name{corr.matrix}
\alias{corr.matrix}
\title{Perform correlations on a matrix, and threshold.}
\usage{
corr.matrix(dat, thresh = NULL, cost = 0.1, exclusions = NULL)
}
\arguments{
  \item{dat}{matrix or data frame of the columns to
  correlate}

  \item{thresh}{absolute correlation value to threshold by}

  \item{cost}{keeps the top X\% of correlations}

  \item{exclusions}{vector of indices (columns) to exclude}
}
\value{
A list with the following components: \item{R}{Pearson
correlation coefficients.} \item{P}{Associated p-values.}
\item{r.thresh}{Binary matrix indicating correlations that
are above a certain threshold.} \item{threshold}{The
threshold used.}
}
\description{
This function does a column-by-column correlation of a
given data frame, and will threshold r-values based on a
given 'cost' value; e.g. 0.1 if you want to keep only the
10\% strongest correlations. It also allows for the
exclusion of a set of columns (i.e. regions or nodes),
given their indices. Also returns the p-values. Essentially
a wrapper for 'rcorr', with some added functionality to
work with the type of data more easily.
}

