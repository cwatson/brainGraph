#' Give vertices in a graph a \emph{lobe} attribute
#'
#' This function will assign vertex attributes \emph{lobe} and \emph{lobe.hemi}
#' for all vertices in a graph, given a specific atlas. It will also add an
#' attribute \emph{circle.layout} for plotting circular graphs.
#'
#' @param g An igraph graph object
#' @param atlas Character string for the atlas name
#' @param atlas.list A list with data for a specific atlas
#'

#' @return An igraph graph object with additional vertex attributes \emph{lobe},
#' \emph{lobe.hemi}, and \emph{circle.layout}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

assign_lobes <- function(g, atlas, atlas.list) {
  lobes <- lobe.hemi <- vector('integer', length=vcount(g))

  lobes[atlas.list$frontal] <- 1
  lobes[atlas.list$parietal] <- 2
  lobes[atlas.list$temporal] <- 3
  lobes[atlas.list$occipital] <- 4
  lobes[atlas.list$insula] <- 5
  lobes[atlas.list$cingulate] <- 6
  
  lobe.hemi[atlas.list$frontal.lh] <- 1
  lobe.hemi[atlas.list$frontal.rh] <- 2
  lobe.hemi[atlas.list$parietal.lh] <- 3
  lobe.hemi[atlas.list$parietal.rh] <- 4
  lobe.hemi[atlas.list$temporal.lh] <- 5
  lobe.hemi[atlas.list$temporal.rh] <- 6
  lobe.hemi[atlas.list$occipital.lh] <- 7
  lobe.hemi[atlas.list$occipital.rh] <- 8
  lobe.hemi[atlas.list$insula.lh] <- 9
  lobe.hemi[atlas.list$insula.rh] <- 10
  lobe.hemi[atlas.list$cingulate.lh] <- 11
  lobe.hemi[atlas.list$cingulate.rh] <- 12

  if (atlas %in% c('dkt', 'dk')) {
    V(g)$circle.layout <- with(atlas.list, c(frontal.lh, insula.lh, cingulate.lh,
                                             temporal.lh, parietal.lh, occipital.lh,
                                             occipital.rh, parietal.rh, temporal.rh,
                                             cingulate.rh, insula.rh, frontal.rh))
    V(g)$hemi[grep('^l.*', rownames(atlas.list$coords))] <- 'L'
    V(g)$hemi[grep('^r.*', rownames(atlas.list$coords))] <- 'R'
    
  } else if (atlas == 'aal90') {
    lobes[atlas.list$limbic] <- 6
    lobes[atlas.list$scgm] <- 7
    V(g)$circle.layout <- with(atlas.list, c(frontal.lh, insula.lh, limbic.lh,
                                             scgm.lh, temporal.lh, parietal.lh,
                                             occipital.lh, occipital.rh,
                                             parietal.rh, temporal.rh, scgm.rh,
                                             limbic.rh, insula.rh, frontal.rh))
    lobe.hemi[atlas.list$limbic.lh] <- 11
    lobe.hemi[atlas.list$limbic.rh] <- 12
    lobe.hemi[atlas.list$scgm.lh] <- 13
    lobe.hemi[atlas.list$scgm.rh] <- 14
  
    V(g)$hemi[grep('.*L$', rownames(atlas.list$coords))] <- 'L'
    V(g)$hemi[grep('.*R$', rownames(atlas.list$coords))] <- 'R'
  
  } else if (atlas %in% c('lpba40', 'hoa112', 'brainsuite')) {
    lobes[atlas.list$scgm] <- 7
    V(g)$circle.layout <- with(atlas.list,
    c(frontal.lh, insula.lh, cingulate.lh, scgm.lh,
    temporal.lh, parietal.lh, occipital.lh,
    occipital.rh, parietal.rh, temporal.rh,
    scgm.rh, cingulate.rh, insula.rh, frontal.rh))
    lobe.hemi[atlas.list$scgm.lh] <- 13
    lobe.hemi[atlas.list$scgm.rh] <- 14
  
    if (atlas == 'lpba40') {
      V(g)$hemi[grep('^l.*', rownames(atlas.list$coords))] <- 'L'
      V(g)$hemi[grep('^r.*', rownames(atlas.list$coords))] <- 'R'
    } else {
      V(g)$hemi[grep('.*L$', rownames(atlas.list$coords))] <- 'L'
      V(g)$hemi[grep('.*R$', rownames(atlas.list$coords))] <- 'R'
    }
  }

  V(g)$lobe <- lobes
  V(g)$lobe.hemi <- lobe.hemi

  g
}
