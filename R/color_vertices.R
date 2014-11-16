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
    lobe.cols <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orange',
                   'lightgreen')

    atlas.list <- eval(parse(text=data(list=atlas)))
    Nv <- nrow(atlas.list$coords)
    lobe.color <- vector(length=Nv)

    lobe.color[atlas.list$frontal] <- lobe.cols[1]
    lobe.color[atlas.list$parietal] <- lobe.cols[2]
    lobe.color[atlas.list$temporal] <- lobe.cols[3]
    lobe.color[atlas.list$occipital] <- lobe.cols[4]
    lobe.color[atlas.list$insula] <- lobe.cols[5]

    if (atlas == 'aal90') {
      lobe.color[atlas.list$limbic] <- lobe.cols[6]
      lobe.color[atlas.list$scgm] <- lobe.cols[7]

    } else if (atlas == 'lpba40' || atlas == 'hoa112' || atlas == 'brainsuite') {
      lobe.color[atlas.list$cingulate] <- lobe.cols[6]
      lobe.color[atlas.list$scgm] <- lobe.cols[7]
    }

    if (atlas == 'lpba40') {
      lobe.color[55:56] <- 'white'
    }

    return(lobe.color)

  } else {
    # Find out how many communities exist that have >= 2 members
    mod.colors <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orangered',
                    'lightgreen', 'lightblue', 'lightyellow')

    mod.sizes <- as.integer(table(comm$membership))
    big.modules <- which(mod.sizes >= 2)
    big.mod.sizes <- mod.sizes[big.modules]
    big.modules <- big.modules[rev(order(big.mod.sizes))]

    mod.colors.comm <- vector(length=length(comm))
    for (i in seq_along(big.modules)) {
      mod.colors.comm[big.modules[i]] <- mod.colors[i]
    }
    mod.colors.comm <- ifelse(mod.colors.comm=='FALSE', 'gray',
                              mod.colors.comm)
 
    return(mod.colors.comm)
  }
}
