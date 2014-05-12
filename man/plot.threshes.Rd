\name{plot.threshes}
\alias{plot.threshes}
\title{Plot vectors for 2 groups along a range of some value.}
\usage{
\method{plot}{threshes}(yvar.1, yvar.2, xlabel = "Cost", ylabel = NULL,
  group1, group2)
}
\arguments{
  \item{yvar.1}{the first group's data}

  \item{yvar.2}{the second group's data}

  \item{xlabel}{the label to place on the x-axis
  (optional)}

  \item{ylabel}{the label to place on the y-axis
  (optional)}

  \item{group1}{Character string identifying first subject
  group}

  \item{group2}{Character string identifying second subject
  group}
}
\description{
This function takes a vector of values from each of two
groups, and plots them over their respective ranges (e.g.
plot the average degree over a range of correlation
thresholds). The default xlabel here is 'Cost'.
}

