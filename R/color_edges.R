#' Color the edges of a graph based on community membership.
#'
#' This function takes the community membership of a given graph, and assigns
#' to the edges a specific color (the same as the vertex membership colors).
#' Currently only works for up to 5 communities (colored: red, green, blue,
#' yellow, magenta). Edges that connect two different communities are colored
#' gray.
#'
#' @param adj.graph The adjacency graph to get its edges colored
#' @param comm The community object returned from community.measures
#' @export
#'

color.edges <- function(adj.graph, comm) {
  mem <- comm$community$membership
  Nc <- max(mem)
  tmp <- list()
  cols <- vector()

  for (i in 1:Nc) {
    if (sum(mem == i) == 1) {
      tmp[[i]] <- 0
    } else {
      tmp[[i]] <- as.vector(E(adj.graph)[which(mem == i) %--% 1:kNumVertices])
      cols[tmp[[i]]] <- comm$vcolors[i]
    }
  }
  dups <- unlist(tmp)[duplicated(unlist(tmp))]
  dups <- dups[dups > 0]

  cols[dups] <- 'gray'
  return(cols)
}
