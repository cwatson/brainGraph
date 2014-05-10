\name{get.resid}
\alias{get.resid}
\title{Run linear models of thickness data to get model residuals.}
\usage{
get.resid(thicknesses, kNumCovars, group1, group2 = NULL)
}
\arguments{
  \item{thicknesses}{data frame of all thickness data,
  including covariates.  First column must be subject ID
  (or something similar)}

  \item{kNumCovars}{Integer specifying number of covariates
  to include in the model}

  \item{group1}{Character string indicating the first
  subject group}

  \item{group2}{Character string indicating the second
  subject group (optional)}
}
\value{
A list with components: \item{models}{The \link{lm} objects
for each brain region.} \item{all}{All residuals for all
subjects and brain regions.} \item{group1}{Residuals for
just group 1.} \item{group2}{Residuals for just group 2.}
}
\description{
This function runs linear models on the thickness data for
each region, in order to adjust for relevant variables
(e.g. age, sex, group, etc.).
}

