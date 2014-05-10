\name{plot.adj}
\alias{plot.adj}
\title{Plots an adjacency matrix with a specific spatial layout.}
\usage{
\method{plot}{adj}(adj.matrix, title = "", subtitle = "",
  vert.lab.col = "black", vert.size = 10, vert.col = "lightblue3",
  vert.lab.cex = 0.75, vert.lab = get.vertex.attribute(adj.matrix, "name"),
  edge.col = "red", edge.width = 2,
  layout = as.matrix(coords.cur[-exclude, 1:2]) \%*\% coord.scale, axes = F,
  rescale = F, xlim = c(-1, 1), ylim = c(-1.5, 1.5))
}
\arguments{
  \item{adj.matrix}{the adjacency matrix to plot}

  \item{vert.label.col}{the color of vertex label text}

  \item{vert.size}{the size of the vertices}

  \item{vert.col}{the color of the vertices}

  \item{vert.lab.cex}{the scaling of vertex label text}

  \item{vert.lab}{the vertex label text}

  \item{edge.col}{the color of the edges}

  \item{edge.width}{the size of the edges}

  \item{layout}{a matrix of the x, y, and z coordinates of
  each vertex}
}
\description{
This function plots an adjacency matrix when the spatial
layout of the nodes is important. It is really just a
wrapper for \code{\link{plot}}. Most of the parameters are
for \code{\link{plot.igraph}}.
}

