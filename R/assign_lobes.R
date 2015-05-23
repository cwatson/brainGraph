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
  atlas <- g$atlas
  vorder <- match(V(g)$name, atlas.dt$name)
  V(g)$lobe <- atlas.dt[vorder, as.numeric(lobe)]
  V(g)$lobe.hemi <- as.numeric(atlas.dt[vorder, interaction(lobe, hemi)])
  V(g)$hemi <- as.character(atlas.dt[vorder, hemi])

  if (atlas %in% c('dkt', 'dk')) {
    third <- 'Cingulate'
    
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
  }

  counts <- atlas.dt[, length(name), by=c('lobe', 'hemi')][order(lobe)]$V1
  V(g)$circle.layout <-
    c(which(V(g)$lobe == 1 & V(g)$hemi == 'L'),
      which(V(g)$lobe == 5 & V(g)$hemi == 'L'),
      which(V(g)$lobe == 6 & V(g)$hemi == 'L'),
      which(V(g)$lobe == 3 & V(g)$hemi == 'L'),
      which(V(g)$lobe == 2 & V(g)$hemi == 'L'),
      which(V(g)$lobe == 4 & V(g)$hemi == 'L'),
      which(V(g)$lobe == 4 & V(g)$hemi == 'R'),
      which(V(g)$lobe == 2 & V(g)$hemi == 'R'),
      which(V(g)$lobe == 3 & V(g)$hemi == 'R'),
      which(V(g)$lobe == 6 & V(g)$hemi == 'R'),
      which(V(g)$lobe == 5 & V(g)$hemi == 'R'),
      which(V(g)$lobe == 1 & V(g)$hemi == 'R'))
#    c(atlas.dt[lobe == 'Frontal' & hemi == 'L', index[vorder]][seq_len(counts[1])],
#      atlas.dt[lobe == 'Insula' & hemi == 'L', index[vorder]][seq_len(counts[9])],
#      atlas.dt[lobe == third & hemi == 'L', index[vorder]][seq_len(counts[11])],
#      atlas.dt[lobe == 'Temporal' & hemi == 'L', index[vorder]][seq_len(counts[5])],
#      atlas.dt[lobe == 'Parietal' & hemi == 'L', index[vorder]][seq_len(counts[3])],
#      atlas.dt[lobe == 'Occipital' & hemi == 'L', index[vorder]][seq_len(counts[7])],
#      atlas.dt[lobe == 'Occipital' & hemi == 'R', index[vorder]][seq_len(counts[8])],
#      atlas.dt[lobe == 'Parietal' & hemi == 'R', index[vorder]][seq_len(counts[4])],
#      atlas.dt[lobe == 'Temporal' & hemi == 'R', index[vorder]][seq_len(counts[6])],
#      atlas.dt[lobe == third & hemi == 'R', index[vorder]][seq_len(counts[12])],
#      atlas.dt[lobe == 'Insula' & hemi == 'R', index[vorder]][seq_len(counts[10])],
#      atlas.dt[lobe == 'Frontal' & hemi == 'R', index[vorder]][seq_len(counts[2])])

  g
}
