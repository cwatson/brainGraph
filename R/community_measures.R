#' Community/module measures.
#'
#' This function splits a graph into communities, and gives vertices for each
#' community a separate color.
#'
#' @param g the adjacency graph
#' @export
#'
#' @return A list with the following components:
#' \item{community}{The \link{communities} object for the present graph.}
#' \item{vcolors}{A list of the vertex colors for the 8 largest communities.}

community.measures <- function(g) {
  comm <- edge.betweenness.community(g)

  # Find out how many communities exist that have >= 2 members
  mod.colors <- c('red', 'green', 'blue', 'yellow', 'magenta', 'lightgreen',
                  'lightblue', 'lightyellow')

  mod.sizes <- vapply(1:length(comm), function(x) sum(comm$membership==x),
                      integer(1))
  big.modules <- which(mod.sizes >= 2)
  big.mod.sizes <- mod.sizes[big.modules]
  big.modules <- big.modules[rev(sort(big.mod.sizes, index.return=T)$ix)]

  mod.colors.comm <- vector(length=length(comm))
  for (i in 1:length(big.modules)) {
    mod.colors.comm[big.modules[i]] <- mod.colors[i]
  }
  mod.colors.comm <- ifelse(mod.colors.comm=='FALSE', 'white',
                              mod.colors.comm)

  return(list(community=comm, vcolors=mod.colors.comm))
}
