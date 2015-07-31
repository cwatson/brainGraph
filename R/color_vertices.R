#' Color graph vertices
#'
#' This function takes an integer vector (representing membership of a community
#' or component) and creates a character vector of colors for each
#' community/module, component, etc. This only assigns a color to groups with at
#' least 2 members; isolated vertices will be colored 'gray'.
#'
#' @param memb An integer vector representing membership of e.g. a community
#'
#' @return A character vector with length equal to the number of communities,
#' lobes, components, etc.

color.vertices <- function(memb) {
  big.modules <- which(as.integer(table(memb)) > 1)

  mod.colors.memb <- rep('gray', length=max(memb))
  mod.colors.memb[big.modules] <- group.cols[big.modules]

  return(mod.colors.memb)
}
