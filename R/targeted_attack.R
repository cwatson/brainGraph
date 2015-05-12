#' Perform a targeted attack of a graph
#'
#' This function performs a "targeted attack" of a graph. It will sort the
#' vertices by either degree or betweenness centrality, and successively remove
#' the top vertices. There is also an option to remove edges successively in
#' terms of edge betweenness. Then it calculates the size of the largest
#' component.
#'
#' @param g The igraph graph object of interest
#' @param type A character string; either 'vertex' or 'edge' removals
#' @param measure A character string; sort by either 'btwn.cent' or 'degree'
#' @export
#'
#' @return A vector representing the percent of maximal component size compared
#' to the graph's original maximal component
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

targeted.attack <- function(g, type=c('vertex', 'edge'),
                            measure=c('btwn.cent', 'degree')) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  if (type == 'vertex') {
    ord <- V(g)$name[order(vertex_attr(g, measure), decreasing=T)]
    max.comp <- vector(length=length(ord)+1)
    max.comp[1] <- g$conn.comp[1, 1]
    for (i in seq_len(vcount(g)-1)) {
      g <- delete.vertices(g, ord[i])
      max.comp[i+1] <- max(components(g)$csize)
    }
  } else {
    if (measure == 'degree') {
      stop('For edge attacks, must choose "btwn.cent"')
    }
    ord <- order(E(g)$btwn, decreasing=T)
    verts <- as_edgelist(g)[ord, ]
    max.comp <- vector(length=length(ord)+1)
    max.comp[1] <- g$conn.comp[1, 1]
    for (i in seq_along(ord)) {
      g <- delete.edges(g, E(g)[V(g)[verts[ord[i], 1]] %--% verts[ord[i], 2]])
      max.comp[i+1] <- max(components(g)$csize)
    }
  }
  return(max.comp / max.comp[1])
}
