#' Give vertices in a graph a \emph{lobe} attribute
#'
#' This function will assign vertex attributes \emph{lobe} and \emph{lobe.hemi}
#' for all vertices in a graph, given a specific atlas. It will also add an
#' attribute \emph{circle.layout} for plotting circular graphs.
#'
#' @param g An igraph graph object
#' @param atlas.dt A list with data for a specific atlas
#'
#' @return An igraph graph object with additional vertex attributes \emph{lobe},
#' \emph{lobe.hemi}, and \emph{circle.layout}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

assign_lobes <- function(g, atlas.dt) {
  lobe <- hemi <- index <- NULL
  vorder <- match(V(g)$name, atlas.dt$name)
  V(g)$lobe <- atlas.dt[vorder, as.numeric(lobe)]
  V(g)$lobe.hemi <- as.numeric(atlas.dt[vorder, interaction(lobe, hemi)])
  V(g)$hemi <- as.character(atlas.dt[vorder, hemi])

  if (atlas %in% c('dkt', 'dk')) {
    third <- 'Cingulate'
    V(g)$circle.layout <- with(atlas.dt, c(frontal.lh, insula.lh, cingulate.lh,
                                             temporal.lh, parietal.lh, occipital.lh,
                                             occipital.rh, parietal.rh, temporal.rh,
                                             cingulate.rh, insula.rh, frontal.rh))
    
  } else if (atlas == 'aal90') {
    V(g)$circle.layout <- with(atlas.dt, c(frontal.lh, insula.lh, limbic.lh,
                                             scgm.lh, temporal.lh, parietal.lh,
                                             occipital.lh, occipital.rh,
                                             parietal.rh, temporal.rh, scgm.rh,
                                             limbic.rh, insula.rh, frontal.rh))
  
  } else if (atlas %in% c('lpba40', 'hoa112', 'brainsuite')) {
    V(g)$circle.layout <- with(atlas.dt, c(frontal.lh, insula.lh, cingulate.lh,
                                             scgm.lh, temporal.lh, parietal.lh,
                                             occipital.lh, occipital.rh,
                                             parietal.rh, temporal.rh, scgm.rh,
                                             cingulate.rh, insula.rh, frontal.rh))
  } else if (atlas == 'destrieux') {
    third <- 'Limbic'
    V(g)$circle.layout <- c(atlas.dt[lobe == 'Frontal' & hemi == 'L', index],
                            atlas.dt[lobe == 'Insula' & hemi == 'L', index],
                            atlas.dt[lobe == third & hemi == 'L', index],
                            atlas.dt[lobe == 'Temporal' & hemi == 'L', index],
                            atlas.dt[lobe == 'Parietal' & hemi == 'L', index],
                            atlas.dt[lobe == 'Occipital' & hemi == 'L', index],
                            atlas.dt[lobe == 'Occipital' & hemi == 'R', index],
                            atlas.dt[lobe == 'Parietal' & hemi == 'R', index],
                            atlas.dt[lobe == 'Temporal' & hemi == 'R', index],
                            atlas.dt[lobe == third & hemi == 'R', index],
                            atlas.dt[lobe == 'Insula' & hemi == 'R', index],
                            atlas.dt[lobe == 'Frontal' & hemi == 'R', index])

  }

  g
}
