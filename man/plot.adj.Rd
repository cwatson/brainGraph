\name{plot.adj}
\alias{plot.adj}
\title{Plots an adjacency matrix with a specific spatial layout.}
\usage{
\method{plot}{adj}(adj.graph, vertex.size = 10, vertex.color = "lightblue3",
  vertex.label.cex = 0.75, edge.color = "red", edge.width = 2,
  layout = as.matrix(coords.cur[-exclude, 1:2]) \%*\% coord.scale,
  rescale = F, ylim = c(-1.5, 1.5), ...)
}
\arguments{
  \item{adj.graph}{The adjacency graph to plot}

  \item{vertex.size}{The size of the vertices}

  \item{vertex.color}{The color of the vertices}

  \item{vertex.label.cex}{The scaling of vertex label text}

  \item{edge.color}{The color of the edges}

  \item{edge.width}{The size of the edges}

  \item{layout}{A matrix of the x, y, (and z, optional)
  coordinates of each vertex}

  \item{rescale}{A logical constant, whether to rescale the
  coordinates}

  \item{ylim}{A vector giving limits for the vertical axis}
}
\description{
This function plots an adjacency matrix when the spatial
layout of the nodes is important. It is really just a
wrapper for \code{\link{plot}}. Most of the parameters are
for \code{\link{plot.igraph}}.
}

