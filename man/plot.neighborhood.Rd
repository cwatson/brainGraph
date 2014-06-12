\name{plot.neighborhood}
\alias{plot.neighborhood}
\title{Plot the neighborhood of a given vertex.}
\usage{
\method{plot}{neighborhood}(g, v, ...)
}
\arguments{
  \item{g}{The graph to be plotted}

  \item{v}{Character string of the vertex name whose
  neighborhood is to be plotted}

  \item{...}{Other parameters (passed to
  \code{\link{plot.adj}})}
}
\description{
This function plots only the neighborhood of a single given
vertex, if you want to see e.g. the neighborhood of the
left precuneus and nothing else.
}

