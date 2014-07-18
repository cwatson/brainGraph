#' Color the edges of a graph based on community membership or major lobe.
#'
#' This function takes the community membership of a given graph, and assigns
#' to the edges a specific color (the same as the vertex membership colors).
#' Currently only works for up to 9 communities (colored: red, green, blue,
#' yellow, magenta, orangered, lightgreen, lightblue, lightyellow). Edges that
#' connect two different communities are colored gray. Also works for the major
#' lobes of the brain (plus insula).
#'
#' @param g The graph to get its edges colored
#' @param comm An integer vector indicating community membership
#' @param lobes An integer vector of the lobe membership for each vertex
#' @param lobe.cols A character vector of the colors each lobe should take
#' @export
#' @return A character vector of colors for each edge in the graph
#'

color.edges <- function(g, comm, lobes=NULL, lobe.cols=NULL) {
  kNumVertices <- vcount(g)
  if (!is.null(lobes)) {
    mem <- lobes
    cols <- lobe.cols
  } else {
    mem <- comm
    cols <- V(g)$color
  }

  kNumComm <- max(mem)
  tmp <- list()
  newcols <- vector()

  for (i in 1:kNumComm) {
    if (sum(mem == i) == 1) {
      tmp[[i]] <- 0
    } else {
      tmp[[i]] <- as.vector(E(g)[which(mem == i) %--% 1:kNumVertices])
      newcols[tmp[[i]]] <- cols[i]
    }
  }
  dups <- unlist(tmp)[duplicated(unlist(tmp))]
  dups <- dups[dups > 0]

  newcols[dups] <- 'gray'

  newcols <- ifelse(is.na(newcols), 'white', newcols)
  
  return(newcols)
}
