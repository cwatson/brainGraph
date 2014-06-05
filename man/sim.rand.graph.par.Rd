\name{sim.rand.graph.par}
\alias{sim.rand.graph.par}
\title{Simulate N random graphs w/ same clustering and degree sequence as the input.}
\usage{
sim.rand.graph.par(g, N)
}
\arguments{
  \item{g}{A graph with the characteristics for simulation
  of random graphs}

  \item{N}{The number of iterations}
}
\value{
A random graph with vertex and graph attributes.
}
\description{
This function will simulate N simple random graphs with the
same clustering and degree sequence as the input. It will
also calculate (and return) avg. path length, transitivity,
mean degree, and modularity. Essentially a wrapper for
\code{\link{sim.rand.graph.clust}} and
\code{\link{set.brainGraph.attributes}}. It uses
\code{\link{foreach}} to speed it up.
}

