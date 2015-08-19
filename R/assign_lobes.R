#' Give vertices in a graph a \emph{lobe} attribute
#'
#' This function will assign vertex attributes \emph{lobe} and \emph{lobe.hemi}
#' for all vertices in a graph, given a specific atlas. It will also add an
#' attribute \emph{circle.layout} for plotting circular graphs.
#'
#' @param g An igraph graph object
#' @param atlas.dt A data table for a specific atlas
#' @param rand A character string indicating whether this function is being run
#' for a random graph or a "graph of interest" (default:FALSE)
#'
#' @return An igraph graph object with additional vertex attributes \emph{lobe},
#' \emph{lobe.hemi}, and \emph{circle.layout}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

assign_lobes <- function(g, atlas.dt, rand=FALSE) {
  lobe <- hemi <- name <- class <- NULL
  atlas <- g$atlas
  vorder <- match(V(g)$name, atlas.dt$name)
  V(g)$lobe <- atlas.dt[vorder, as.numeric(lobe)]
  V(g)$lobe.hemi <- as.numeric(atlas.dt[vorder, interaction(lobe, hemi)])
  V(g)$hemi <- as.character(atlas.dt[vorder, hemi])

  if (atlas == 'destrieux') V(g)$class <- atlas.dt[vorder, as.numeric(class)]

  if (!isTRUE(rand)) {
    if (atlas %in% c('dkt', 'dk', 'destrieux')) {
      counts <- atlas.dt[order(lobe), .N, by=.(lobe, hemi)]$N
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

    } else if (atlas %in% c('aal90', 'lpba40', 'hoa112', 'brainsuite',
                            'dk.scgm', 'dkt.scgm')) {
      counts <- atlas.dt[order(lobe), .N, by=.(lobe, hemi)]$N
      V(g)$circle.layout <-
        c(which(V(g)$lobe == 1 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 5 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 6 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 7 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 3 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 2 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 4 & V(g)$hemi == 'L'),
          which(V(g)$lobe == 4 & V(g)$hemi == 'R'),
          which(V(g)$lobe == 2 & V(g)$hemi == 'R'),
          which(V(g)$lobe == 3 & V(g)$hemi == 'R'),
          which(V(g)$lobe == 7 & V(g)$hemi == 'R'),
          which(V(g)$lobe == 6 & V(g)$hemi == 'R'),
          which(V(g)$lobe == 5 & V(g)$hemi == 'R'),
          which(V(g)$lobe == 1 & V(g)$hemi == 'R'))
    }
  }

  g
}
