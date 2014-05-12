\name{sim.rand.graph.par}
\alias{sim.rand.graph.par}
\title{Simulate N random graphs w/ same clustering and degree sequence as the input.}
\usage{
sim.rand.graph.par(deg.dist, N, cl)
}
\arguments{
  \item{deg.dist}{A vector of the degree distribution the
  random graphs should have}

  \item{N}{The number of iterations}

  \item{cl}{The clustering (transitivity) the random graphs
  should have}
}
\value{
A list with the following components: \item{Lp.rand}{A
vector of the characteristic path lengths for the N random
graphs that were created.} \item{Cp.rand}{A vector of the
clustering coefficients for the N random graphs that were
created.} \item{g.eff}{A vector of the global efficiencies
for the N random graphs that were created.} \item{l.eff}{A
vector of the mean local efficiencies for the N random
graphs that were created.} \item{rich}{A vector of the rich
club coefficients for the N random graphs that were
created.} \item{mod}{A vector of the modularity for the N
random graphs.} \item{Ek}{The number of edges for each
graph.} \item{triangles}{The number of triangles for each
graph.}
}
\description{
This function will simulate N simple random graphs with the
same clustering and degree sequence as the input. It will
also calculate (and return) avg. path length, transitivity,
mean degree, and modularity. Essentially a wrapper for
\code{\link{sim.rand.graph.clust}}. It uses
\code{\link{foreach}} to speed it up.
}

