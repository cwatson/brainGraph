#' Color graph edges
#'
#' This function takes the community membership of a given graph, and assigns
#' to the edges a specific color (the same as the vertex membership colors).
#' Edges that connect vertices of two different groups are colored gray. Also
#' works for the major lobes of the brain, plus insula, subcortical gray matter,
#' cingulate, limbic lobe (if included in the specific brain atlas).
#'
#' @param g The graph to get its edges colored
#' @param memb An integer vector indicating vertex group membership
#'
#' @return A character vector of colors for each edge in the graph

color.edges <- function(g, memb) {
  big.modules <- as.integer(names(which(table(memb) > 1)))

  newcols <- rep('gray50', length=ecount(g))
  tmp <- vector('list', length=max(big.modules))
  for (i in big.modules) {
    x <- which(memb == i)
    tmp[[i]] <- as.vector(E(g)[x %--% x])
  }

  for (i in big.modules) {
    if (!is.null(tmp[[i]])) {
      newcols[tmp[[i]]] <- group.cols[i]
    }
  }

  return(newcols)
}
