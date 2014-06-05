\name{set.brainGraph.attributes}
\alias{set.brainGraph.attributes}
\title{Set a number of graph and vertex attributes useful in MRI analyses}
\usage{
set.brainGraph.attributes(g, coords = NULL, rand = FALSE)
}
\arguments{
  \item{g}{A graph}

  \item{coords}{A matrix of the spatial coordinates}

  \item{rand}{Logical indicating whether the graph is
  random or not}
}
\value{
g A copy of the same graph, with the attributes added.
}
\description{
This function will set a number of graph and vertex
attributes of a given igraph object. It is required that a
community structure has already been determined, and
\code{\link{part.coeff}} and
\code{\link{within.module.deg.z.score}} have been
calculated.
}

