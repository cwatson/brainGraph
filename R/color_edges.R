#' Color the edges of a graph based on community membership.
#'
#' This function takes the community membership of a given graph, and assigns
#' to the edges a specific color (the same as the vertex membership colors).
#' Currently only works for up to 9 communities (colored: red, green, blue,
#' yellow, magenta, orangered, lightgreen, lightblue, lightyellow). Edges that
#' connect two different communities are colored gray.
#'
#' @param g The graph to get its edges colored
#' @param comm The community object returned from community.measures
#' @export
#' @return A character vector of colors for each edge in the graph
#'

color.edges <- function(g, comm) {
  kNumVertices <- vcount(g)
  mem <- comm$community$membership
  kNumComm <- max(mem)
  tmp <- list()
  cols <- vector()

  for (i in 1:kNumComm) {
    if (sum(mem == i) == 1) {
      tmp[[i]] <- 0
    } else {
      tmp[[i]] <- as.vector(E(g)[which(mem == i) %--% 1:kNumVertices])
      cols[tmp[[i]]] <- comm$vcolors[i]
    }
  }
  dups <- unlist(tmp)[duplicated(unlist(tmp))]
  dups <- dups[dups > 0]

  cols[dups] <- 'gray'

  cols <- ifelse(is.na(cols), 'white', cols)
  return(cols)
}
