#' Select edges for re-wiring.
#'
#' This function selects edges to be re-wired when simulating random graphs
#' while controlling for \emph{clustering}. It is based on the algorithm by
#' Bansal et al. (2009), BMC Bioinformatics.
#'
#' @param A Numeric (adjacency) matrix
#' @param degs.large Integer vector of vertex numbers with degree greater than
#'   one
#'
#' @return A data frame with four elements; two edges will be removed and two
#'   edges will be added between the four vertices.
#'
#' @family Null graph functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Bansal S., Khandelwal S., Meyers L.A. (2009) \emph{Exploring
#'   biological network structure with clustered random networks}. BMC
#'   Bioinformatics, 10:405-421.

choose.edges <- function(A, degs.large) {

  # Uniformly select random node with degree > 1
  #=============================================================================
  get.node <- function(A, degs.large) {
    repeat {
      x <- sample(degs.large, 1)
      neighb <- intersect(which(A[x, ] == 1), degs.large)
      if (length(neighb) >= 2) return(list(x, neighb))
    }
  }
  # Uniformly select 2 random neighbors with degree > 1
  #=============================================================================
  get.neighbors <- function(nbrhood) {
    y <- sample(nbrhood, 2)
    return(data.frame(y1=y[1], y2=y[2]))
  }
  #=============================================================================
  # Uniformly select a random neighbor from the y's with degree > 1
  #=============================================================================
  get.neighbors.z1 <- function(A, node, y1) {
    y1.neighb <- which(A[y1, ] == 1)
    choices <- setdiff(y1.neighb, node)
    if (length(choices) <= 1) {
      return(choices)
    } else {
      return(sample(choices, 1))
    }
  }
  #=============================================================================
  get.neighbors.z2 <- function(A, node, degs.large, y2, z1) {
    y2.neighb <- which(A[y2, ] == 1)
    repeat {
      choices <- setdiff(y2.neighb, c(node, z1))
      if (length(choices) == 1) {
        return(choices)
      } else if (length(choices) > 1) {
        return(sample(choices, 1))
      }
      tmp <- get.node(A, degs.large)
      node <- tmp[[1]]
      neighb <- tmp[[2]]
      n <- get.neighbors(neighb)
      y1 <- n$y1
      y2 <- n$y2

      z1 <- get.neighbors.z1(A, node, y1)
      get.neighbors.z2(A, node, degs.large, y2, z1)
    }
  }
  #=============================================================================
  #=============================================================================
  tmp <- get.node(A, degs.large)
  x <- tmp[[1]]
  neighb <- tmp[[2]]

  n <- get.neighbors(neighb)
  y1 <- n$y1
  y2 <- n$y2
  z1 <- get.neighbors.z1(A, x, y1)
  z2 <- get.neighbors.z2(A, x, degs.large, y2, z1)

  return(data.frame(y1, y2, z1, z2))
}
