#' Calculate the rich core of a graph
#'
#' This function finds the boundary of the rich core of a graph, based on the
#' decreasing order of vertex degree. It also calculates the degree that
#' corresponds that that rank, and the core size relative to the total number of
#' vertices in the graph.
#'
#' @param g The graph of interest
#' @export
#'
#' @return A data frame with the following components:
#' \item{density}{The density of the graph.}
#' \item{rank}{The rank of the boundary for the rich core.}
#' \item{k.r}{The degree of the vertex at the boundary.}
#' \item{core.size}{The size of the core relative to the graph size.}
#'
#' @seealso \code{\link{rich.club.coeff}, \link{rich.club.norm}}
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Ma A \& Mondragon R.J. (2015) \emph{Rich-cores in networks}. PLoS
#' One, 10(3): e0119678. doi:10.1371/journal.pone.0119678

rich.core <- function(g) {
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  if ('density' %in% graph_attr_names(g)) {
    dens <- g$density
  } else {
    dens <- graph.density(g)
  }
  Nv <- vcount(g)

  vorder <- order(degs, decreasing=TRUE)
  kplus <- sapply(seq_len(Nv), function(x)
                  length(E(g)[vorder[x] %--% which(degs > degs[vorder[x]])]))

  r <- max(which(kplus == max(kplus)))
  k.r <- degs[vorder][r]
  core.size <- r / Nv
  return(data.frame(density=dens, rank=r, k.r=k.r, core.size=core.size))
}
