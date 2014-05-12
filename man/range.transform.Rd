\name{range.transform}
\alias{range.transform}
\title{Transform a vector to have a different range.}
\usage{
\method{range}{transform}(x, min.val = 0, max.val = 1)
}
\arguments{
  \item{x}{the vector to transform}

  \item{min.val}{the minimum value of the new range}

  \item{max.val}{the maximum value of the new range}
}
\value{
A vector of the transformed input.
}
\description{
This function takes a vector and transforms it to have a
new range, given the input, or the default values of [0,
1].
}

