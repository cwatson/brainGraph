#' Color graph vertices
#'
#' This function takes an integer vector (representing membership of a community
#' or component) and creates a character vector of colors for each
#' community/module, component, etc.
#'
#' @param memb An integer vector representing membership of e.g. a community
#' @param cols A character vector of the colors each vertex group should take
#'
#' @return A character vector with length equal to the number of communities,
#' lobes, components, etc.

color.vertices <- function(memb, cols) {

  # Find out how many communities exist that have >= 2 members
  mod.sizes <- as.integer(table(memb))
  big.modules <- which(mod.sizes >= 2)
  big.mod.sizes <- mod.sizes[big.modules]
  big.modules <- big.modules[rev(order(big.mod.sizes))]

  mod.colors.memb <- vector('character', length=max(memb))
  for (i in seq_along(big.modules)) {
    mod.colors.memb[big.modules[i]] <- cols[i]
  }
  mod.colors.memb <- ifelse(mod.colors.memb=='FALSE', 'gray',
                            mod.colors.memb)

  return(mod.colors.memb)
}
