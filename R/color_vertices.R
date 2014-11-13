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

    } else if (atlas == 'aal90') {
      Nv <- nrow(aal90$coords)
      lobe.color <- vector(length=Nv)
      lobe.color[aal90$frontal] <- lobe.cols[1]
      lobe.color[aal90$parietal] <- lobe.cols[2]
      lobe.color[aal90$temporal] <- lobe.cols[3]
      lobe.color[aal90$occipital] <- lobe.cols[4]
      lobe.color[aal90$insula] <- lobe.cols[5]
      lobe.color[aal90$limbic] <- lobe.cols[6]
      lobe.color[aal90$scgm] <- lobe.cols[7]

    } else if (atlas == 'lpba40') {
      Nv <- nrow(lpba40$coords)
      lobe.color <- vector(length=Nv)
      lobe.color[lpba40$frontal] <- lobe.cols[1]
      lobe.color[lpba40$parietal] <- lobe.cols[2]
      lobe.color[lpba40$temporal] <- lobe.cols[3]
      lobe.color[lpba40$occipital] <- lobe.cols[4]
      lobe.color[lpba40$insula] <- lobe.cols[5]
      lobe.color[lpba40$cingulate] <- lobe.cols[6]
      lobe.color[lpba40$scgm] <- lobe.cols[7]
      lobe.color[55:56] <- 'white'

    } else if (atlas == 'hoa112') {
      Nv <- nrow(hoa112$coords)
      lobe.color <- vector(length=Nv)
      lobe.color[hoa112$frontal] <- lobe.cols[1]
      lobe.color[hoa112$parietal] <- lobe.cols[2]
      lobe.color[hoa112$temporal] <- lobe.cols[3]
      lobe.color[hoa112$occipital] <- lobe.cols[4]
      lobe.color[hoa112$insula] <- lobe.cols[5]
      lobe.color[hoa112$cingulate] <- lobe.cols[6]
      lobe.color[hoa112$scgm] <- lobe.cols[7]

    } else if (atlas == 'brainsuite') {
      Nv <- nrow(brainsuite$coords)
      lobe.color <- vector(length=Nv)
      lobe.color[brainsuite$frontal] <- lobe.cols[1]
      lobe.color[brainsuite$parietal] <- lobe.cols[2]
      lobe.color[brainsuite$temporal] <- lobe.cols[3]
      lobe.color[brainsuite$occipital] <- lobe.cols[4]
      lobe.color[brainsuite$insula] <- lobe.cols[5]
      lobe.color[brainsuite$cingulate] <- lobe.cols[6]
      lobe.color[brainsuite$scgm] <- lobe.cols[7]
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
    for (i in 1:length(big.modules)) {
      mod.colors.comm[big.modules[i]] <- mod.colors[i]
    }
    mod.colors.comm <- ifelse(mod.colors.comm=='FALSE', 'gray',
                              mod.colors.comm)
 
    return(mod.colors.comm)
  }
}
