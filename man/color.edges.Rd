\name{color.edges}
\alias{color.edges}
\title{Color the edges of a graph based on community membership.}
\usage{
color.edges(g, comm)
}
\arguments{
  \item{g}{The adjacency graph to get its edges colored}

  \item{comm}{The community object returned from
  community.measures}
}
\description{
This function takes the community membership of a given
graph, and assigns to the edges a specific color (the same
as the vertex membership colors). Currently only works for
up to 5 communities (colored: red, green, blue, yellow,
magenta). Edges that connect two different communities are
colored gray.
}

