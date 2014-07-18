#' Give graph vertices a color attribute.
#'
#' This function takes a \link{communities} object (or a character string
#' indicating the brain atlas used (e.g. for the major lobes of the brain), and
#' creates a character vector of colors for each vertex.
#'
#' @param comm A \link{communities} object
#' @param atlas A character string
#' @export
#'
#' @return A character vector of length \code{vcount(g)}

color.vertices <- function(comm, atlas=NULL) {
  if (!is.null(atlas)) {
    lobe.cols <- c('red', 'green', 'blue', 'magenta', 'yellow')
    if (atlas == 'dkt') { 
      Nv <- nrow(dkt$coords)
      lobe.color <- vector(length=Nv)
      lobe.color[dkt$frontal] <- lobe.cols[1]
      lobe.color[dkt$parietal] <- lobe.cols[2]
      lobe.color[dkt$temporal] <- lobe.cols[3]
      lobe.color[dkt$occipital] <- lobe.cols[4]
      lobe.color[dkt$insula] <- lobe.cols[5]

    } else if (atlas == 'dk') { 
      Nv <- nrow(dk$coords)
      lobe.color <- vector(length=Nv)
      lobe.color[dk$frontal] <- lobe.cols[1]
      lobe.color[dk$parietal] <- lobe.cols[2]
      lobe.color[dk$temporal] <- lobe.cols[3]
      lobe.color[dk$occipital] <- lobe.cols[4]
      lobe.color[dk$insula] <- lobe.cols[5]
    }
    return(lobe.color)

  } else {
    # Find out how many communities exist that have >= 2 members
    mod.colors <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orangered',
                    'lightgreen', 'lightblue', 'lightyellow')

    mod.sizes <- as.integer(table(comm$membership))
    big.modules <- which(mod.sizes >= 2)
    big.mod.sizes <- mod.sizes[big.modules]
    big.modules <- big.modules[rev(sort(big.mod.sizes, index.return=T)$ix)]

    mod.colors.comm <- vector(length=length(comm))
    for (i in 1:length(big.modules)) {
      mod.colors.comm[big.modules[i]] <- mod.colors[i]
    }
    mod.colors.comm <- ifelse(mod.colors.comm=='FALSE', 'white',
                              mod.colors.comm)
 
    return(mod.colors.comm)
  }
}
