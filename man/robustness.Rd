% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robustness.R
\name{robustness}
\alias{robustness}
\title{Analysis of network robustness}
\usage{
robustness(g, type = c("vertex", "edge"), measure = c("btwn.cent",
  "degree", "random"), N = 1000)
}
\arguments{
\item{g}{An \code{igraph} graph object}

\item{type}{Character string; either \code{'vertex'} or \code{'edge'}
removals. Default: \code{vertex}}

\item{measure}{Character string; sort by either \code{'btwn.cent'} or
\code{'degree'}, or choose \code{'random'}. Default: \code{'btwn.cent'}}

\item{N}{Integer; the number of iterations if \code{'random'} is chosen.
Default: \code{1e3}}
}
\value{
Data table with elements:
  \item{type}{Character string describing the type of analysis performed}
  \item{measure}{The input argument \code{measure}}
  \item{comp.size}{The size of the largest component after edge/vertex
    removal}
  \item{comp.pct}{Numeric vector of the ratio of maximal component size after
    each removal to the observed graph's maximal component size}
  \item{removed.pct}{Numeric vector of the ratio of vertices/edges removed}
  \item{Group}{Character string indicating the subject group, if applicable}
}
\description{
This function performs a \dQuote{targeted attack} of a graph or a
\dQuote{random failure} analysis, calculating the size of the largest
component after edge or vertex removal.
}
\details{
In a targeted attack, it will sort the vertices by either degree or
betweenness centrality (or sort edges by betweenness), and successively
remove the top vertices/edges. Then it calculates the size of the largest
component.

In a random failure analysis, vertices/edges are removed in a random order.
}
\references{
Albert, R. and Jeong, H. and Barabasi, A. (2000) Error and attack
  tolerance of complex networks. \emph{Nature}, \bold{406}, 378--381.
  \doi{10.1038/35019019}
}
\author{
Christopher G. Watson, \email{cgwatson@bu.edu}
}
