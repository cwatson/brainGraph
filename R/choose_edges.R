#' Select edges for re-wiring.
#'
#' This function selects edges to be re-wired when simulating random graphs. It
#' is based on the algorithm by Bansal et al. (2009), BMC Bioinformatics.
#'
#' @param g The random graph that has been generated
#'
#' @return A data frame with four elements; two edges will be removed and two
#' edges will be added between the four vertices.
#'
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#' biological network structure with clustered random networks}. BMC
#' Bioinformatics, 10:405-421.

choose.edges <- function(g) {
  degs <- degree(g)
  degs.large <- as.integer(which(degs > 1))
   
  #=============================================================================
  # Uniformly select random node with degree > 1
  #=============================================================================
  get.node <- function(graph, degrees.large) {
    repeat {
      x <- sample(degrees.large, 1)
      neighb <- intersect(neighbors(graph, x), degrees.large)
      if (length(neighb) >= 2) return(list(x, neighb))
    }
  }
  #=============================================================================
  # Uniformly select 2 random neighbors with degree > 1
  #=============================================================================
  get.neighbors <- function(nbrhood) {
    y <- sample(nbrhood, 2)
    return(data.frame(y1=y[1], y2=y[2]))
  }
  #=============================================================================
  # Uniformly select a random neighbor from the y's with degree > 1
  #=============================================================================
  get.neighbors.z1 <- function(graph, node, degrees, y1) {
    y1.neighb <- as.integer(neighbors(graph, y1))
    choices <- setdiff(y1.neighb, node)
    if (length(choices) <= 1) {
      return(choices)
    } else {
      return(sample(choices, 1))
    }
  }
  #=============================================================================
  get.neighbors.z2 <- function(graph, node, degrees, degrees.large, y2, z1) {
    #iter <- 0
    y2.neighb <- as.integer(neighbors(graph, y2))
    repeat {
      choices <- setdiff(y2.neighb, c(node, z1))
      if (length(choices) == 1) {
        return(choices)
      } else if (length(choices) > 1) {
        return(sample(choices, 1))
      }
      #iter <- iter + 1
      #if (iter >= length(y2.neighb)) {
        tmp <- get.node(graph, degrees.large)
        node <- tmp[[1]]
        neighb <- tmp[[2]]
        n <- get.neighbors(neighb)
        y1 <- n$y1
        y2 <- n$y2

        z1 <- get.neighbors.z1(graph, node, degrees, y1)
        get.neighbors.z2(graph, node, degrees, degrees.large, y2, z1)
      #}
    }
  }
  #=============================================================================
  #=============================================================================
  tmp <- get.node(g, degs.large)
  x <- tmp[[1]]
  neighb <- tmp[[2]]

  n <- get.neighbors(neighb)
  y1 <- n$y1
  y2 <- n$y2
  z1 <- get.neighbors.z1(g, x, degs, y1)
  z2 <- get.neighbors.z2(g, x, degs, degs.large, y2, z1)

  return(data.frame(x, y1, y2, z1, z2))
}
