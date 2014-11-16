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

  # Color edges based on "lobe"
  if (!is.null(lobes)) {
    mem <- lobes
    kNumComm <- max(mem)
    comm.order <- 1:kNumComm
    cols <- lobe.cols
  } else {
  # Color edges based on community membership
    mem <- comm
    kNumComm <- max(mem)
    comm.order <- rev(order(as.integer(table(mem))))
    cols <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orangered',
              'lightgreen', 'lightblue', 'lightyellow')
  }

  tmp <- vector('list', length=kNumComm)
  newcols <- vector('character', length=ecount(g))

  for (i in seq_len(kNumComm)) {
    if (sum(mem == comm.order[i]) == 1) {
      tmp[[i]] <- 0
    } else {
      tmp[[i]] <- as.vector(E(g)[which(mem == comm.order[i]) %--% which(mem == comm.order[i])])
      newcols[tmp[[i]]] <- cols[i]
    }
  }
  dups <- unlist(tmp)[duplicated(unlist(tmp))]
  dups <- dups[dups > 0]

  newcols[dups] <- 'gray50'

  newcols <- ifelse(newcols=='', 'gray50', newcols)
  
  return(newcols)
}
