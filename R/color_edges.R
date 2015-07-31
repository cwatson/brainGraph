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
  kNumComm <- max(memb)
  comm.order <- seq_len(kNumComm)
  cols <- group.cols

  tmp <- vector('list', length=kNumComm)
  newcols <- rep('gray50', length=ecount(g))

  sums <- vapply(seq_len(kNumComm), function(x) sum(memb == comm.order[x]), integer(1))
  tmp[which(sums == 1)] <- 0
  for (i in which(sums > 1)) {
      matches <- which(memb == comm.order[i])
      tmp[[i]] <- as.vector(E(g)[matches %--% matches])
      newcols[tmp[[i]]] <- cols[i]
  }

  return(newcols)
}
