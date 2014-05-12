\name{rich.club.coeff}
\alias{rich.club.coeff}
\title{Calculate the rich club of a graph.}
\usage{
rich.club.coeff(graph, k)
}
\arguments{
  \item{graph}{the adjacency graph of interest}

  \item{k}{the minimum degree for including a node}
}
\value{
A list with the following components: \item{coeff}{The rich
club coefficient, phi.} \item{graph}{A subgraph containing
only the rich club nodes.} \item{Nk}{The number of vertices
in the rich club graph.} \item{Ek}{The number of edges in
the rich club graph.}
}
\description{
This function calculates the rich club of a graph, both the
coefficient phi and the nodes that make up this subgraph.
}

