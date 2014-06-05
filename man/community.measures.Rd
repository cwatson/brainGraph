\name{community.measures}
\alias{community.measures}
\title{Community/module measures.}
\usage{
community.measures(g)
}
\arguments{
  \item{g}{the adjacency graph}
}
\value{
A list with the following components: \item{community}{The
\link{communities} object for the present graph.}
\item{vcolors}{A list of the vertex colors for the 8
largest communities.}
}
\description{
This function splits a graph into communities, and gives
vertices for each community a separate color. Any "small"
communities (defined by me as having degree < 2) are
colored white.
}

