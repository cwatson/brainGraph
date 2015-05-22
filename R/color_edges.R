#' Color graph edges
#'
#' This function takes the community membership of a given graph, and assigns
#' to the edges a specific color (the same as the vertex membership colors).
#' Currently only works for up to 9 communities (colored: red, green, blue,
#' yellow, magenta, orange, lightgreen, lightblue, lightyellow). Edges that
#' connect two different communities are colored gray. Also works for the major
#' lobes of the brain (plus insula).
#'
#' @param g The graph to get its edges colored
#' @param comm An integer vector indicating community membership
#' @param lobes An integer vector of the lobe membership for each vertex
#' @param cols A character vector of the colors each vertex group should take
#'
#' @return A character vector of colors for each edge in the graph

color.edges <- function(g, comm, lobes=NULL, cols=NULL) {
  # Color edges based on "lobe"
  if (!is.null(lobes)) {
    mem <- lobes
    kNumComm <- max(mem)
    comm.order <- seq_len(kNumComm)
  } else {
  # Color edges based on community membership
    mem <- comm
    kNumComm <- max(mem)
    comm.order <- rev(order(as.integer(table(mem))))
  }

  tmp <- vector('list', length=kNumComm)
  newcols <- vector('character', length=ecount(g))

  sums <- sapply(seq_len(kNumComm), function(x) sum(mem == comm.order[x]))
  tmp[which(sums == 1)] <- 0
  for (i in which(sums > 1)) {
      matches <- which(mem == comm.order[i])
      tmp[[i]] <- as.vector(E(g)[matches %--% matches])
      newcols[tmp[[i]]] <- cols[i]
  }

  newcols <- ifelse(newcols=='', 'gray50', newcols)
}
