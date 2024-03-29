% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_brainnet.R
\name{write_brainnet}
\alias{write_brainnet}
\title{Write files to be used for visualization with BrainNet Viewer}
\usage{
write_brainnet(g, vcolor = "none", vsize = "constant",
  edge.wt = NULL, file.prefix = "")
}
\arguments{
\item{g}{The \code{igraph} graph object of interest}

\item{vcolor}{Character string indicating how to color the vertices (default:
\code{'none'})}

\item{vsize}{Character string indicating what size the vertices should be;
can be any vertex-level attribute (default: \code{'constant'})}

\item{edge.wt}{Character string indicating the edge attribute to use to
return a weighted adjacency matrix (default: \code{NULL})}

\item{file.prefix}{Character string for the basename of the \code{.node} and
\code{.edge} files that are written}
}
\description{
Write the \code{.node} and \code{.edge} files necessary for visualization
with the BrainNet Viewer software.
}
\details{
For the \code{.node} file, there are 6 columns:
\itemize{
  \item \emph{Columns 1-3}: Vertex x-, y-, and z-coordinates
  \item \emph{Column 4}: Vertex color
  \item \emph{Column 5}: Vertex size
  \item \emph{Column 6}: Vertex label
}
The \code{.edge} file is the graph's associated adjacency matrix; a weighted
adjacency matrix can be returned by using the \code{edge.wt} argument.
}
\examples{
\dontrun{
write_brainnet(g, vcolor='community', vsize='degree', edge.wt='t.stat')
}
}
\references{
Xia, M. and Wang, J. and He, Y. (2013). BrainNet Viewer: a
  network visualization tool for human brain connectomics. \emph{PLoS One},
  \bold{8(7)}, e68910. \doi{10.1371/journal.pone.0068910}
}
\author{
Christopher G. Watson, \email{cgwatson@bu.edu}
}
