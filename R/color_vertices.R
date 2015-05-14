#' Color graph vertices
#'
#' This function takes an integer vector (representing membership of a community
#' or component), or a character string indicating the brain atlas used (e.g. for
#' the major lobes of the brain), and creates a character vector of colors for
#' each community/module, component, etc.
#'
#' @param memb An integer vector representing membership of e.g. a community
#' @param atlas A character string
#'
#' @return A character vector with length equal to the number of communities,
#' lobes, components, etc.

color.vertices <- function(memb, atlas=NULL) {
  cols <- c('red', 'green', 'blue', 'magenta', 'yellow', 'orange',
            'lightgreen', 'lightblue', 'lightyellow')

  if (!is.null(atlas)) {
    atlas.list <- eval(parse(text=atlas))
    Nv <- nrow(atlas.list$coords)
    lobe.color <- vector(length=Nv)

    lobe.color[atlas.list$frontal] <- cols[1]
    lobe.color[atlas.list$parietal] <- cols[2]
    lobe.color[atlas.list$temporal] <- cols[3]
    lobe.color[atlas.list$occipital] <- cols[4]
    lobe.color[atlas.list$insula] <- cols[5]
    lobe.color[atlas.list$cingulate] <- cols[6]

    if (atlas == 'aal90') {
      lobe.color[atlas.list$limbic] <- cols[6]
      lobe.color[atlas.list$scgm] <- cols[7]

    } else if (atlas %in% c('lpba40', 'hoa112', 'brainsuite')) {
      lobe.color[atlas.list$scgm] <- cols[7]

    } else if (atlas == 'destrieux') {
      lobe.color[atlas.list$limbic] <- cols[6]
    }

    if (atlas == 'lpba40') {
      lobe.color[55:56] <- 'white'
    }

    return(lobe.color)

  } else {
    # Find out how many communities exist that have >= 2 members
    mod.sizes <- as.integer(table(memb))
    big.modules <- which(mod.sizes >= 2)
    big.mod.sizes <- mod.sizes[big.modules]
    big.modules <- big.modules[rev(order(big.mod.sizes))]

    mod.colors.memb <- vector(length=max(memb))
    for (i in seq_along(big.modules)) {
      mod.colors.memb[big.modules[i]] <- cols[i]
    }
    mod.colors.memb <- ifelse(mod.colors.memb=='FALSE', 'gray',
                              mod.colors.memb)
 
    return(mod.colors.memb)
  }
}
