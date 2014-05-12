\name{part.coeff}
\alias{part.coeff}
\title{Calculate the participation coefficient of each node}
\usage{
part.coeff(adj.mat, adj.graph, community.mem)
}
\arguments{
  \item{adj.mat}{the adjacency matrix}

  \item{adj.graph}{the adjacency graph}

  \item{community.mem}{the community membership indices of
  each node}
}
\value{
A vector of the participation coeff's for each node of the
graph.
}
\description{
This function calculates the participation coefficient of
each node in a graph. The graph must be split into modules
beforehand. The coefficient should not exceed 1 - (1 / #
modules).
}
\references{
Guimera, R. and Amaral, L.A.N. (2005) Cartography of
complex networks: modules and universal roles, Journal of
Statistical Mechanics: Theory and Experiment, 02, P02001.
}

