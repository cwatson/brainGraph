#' Give vertices in a graph a \emph{lobe} attribute.
#'
#' This function will assign vertex attributes \emph{lobe} and \emph{lobe.hemi}
#' for all vertices in a graph, given a specific atlas. It will also add an
#' attribute \emph{circle.layout} for plotting circular graphs.
#'
#' The input graph \code{g} \emph{must} have a graph attribute named \code{atlas}.
#'
#' @param g An \emph{igraph} graph object.
#' @param rand A character string indicating whether this function is being run
#' for a random graph or a "graph of interest" (default: FALSE).
#'
#' @return An \emph{igraph} graph object with additional vertex attributes:
#'   \item{lobe}{Character string indicating the lobe}
#'   \item{lobe.hemi}{Integer vector indicating the lobe and hemisphere}
#'   \item{circle.layout}{Integer vector for ordering the vertices for circle
#'     plots}
#'   \item{x, y, z, x.mni, y.mni, z.mni}{Spatial coordinates}
#'   \item{color.lobe}{Colors based on \emph{lobe}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

assign_lobes <- function(g, rand=FALSE) {

  stopifnot(is_igraph(g))
  stopifnot('atlas' %in% graph_attr_names(g))
  lobe <- hemi <- name <- class <- x <- y <- z <- x.mni <- y.mni <- z.mni <- NULL

  # Check that vertex names match the atlas names
  atlas.dt <- eval(parse(text=g$atlas))
  nonmatches <- !V(g)$name %in% atlas.dt[, name]
  if (any(nonmatches)) {
    stop(paste('Check the following vertex names: ',
               paste(V(g)$name[nonmatches], collapse=' ')))
  }
  vorder <- match(V(g)$name, atlas.dt$name)

  V(g)$lobe <- atlas.dt[vorder, as.numeric(lobe)]
  V(g)$lobe.hemi <- as.numeric(atlas.dt[vorder, interaction(lobe, hemi)])
  V(g)$hemi <- as.character(atlas.dt[vorder, hemi])

  if (g$atlas == 'destrieux') V(g)$class <- atlas.dt[vorder, as.numeric(class)]

  if (!isTRUE(rand)) {
    # Add spatial coordinates for plotting over a brain slice
    V(g)$x <- atlas.dt[vorder, x]
    V(g)$y <- atlas.dt[vorder, y]
    V(g)$z <- atlas.dt[vorder, z]
    V(g)$x.mni <- atlas.dt[vorder, x.mni]
    V(g)$y.mni <- atlas.dt[vorder, y.mni]
    V(g)$z.mni <- atlas.dt[vorder, z.mni]
    V(g)$color.lobe <- group.cols[V(g)$lobe]
    E(g)$color.lobe <- color.edges(g, V(g)$lobe)

    counts <- atlas.dt[order(lobe), .N, by=list(lobe, hemi)]$N
    if (g$atlas %in% c('dkt', 'dk', 'destrieux')) {
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

    } else if (g$atlas %in% c('aal90', 'aal2.94', 'aal116', 'aal2.120', 'lpba40',
                              'hoa112', 'brainsuite', 'dk.scgm', 'dkt.scgm')) {
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
      if (g$atlas %in% c('aal116', 'aal2.120')) {
        mid1 <- sum(counts[seq(1, 14, 2)])
        mid2 <- sum(counts[seq(1, 15, 2)])
        V(g)$circle.layout <-
          c(V(g)$circle.layout[1:mid1],
            which(V(g)$lobe == 8 & V(g)$hemi == 'L'),
            V(g)$circle.layout[(mid2+1):(mid2+mid1)],
            which(V(g)$lobe == 8 & V(g)$hemi == 'R'))
      }
    }
  }

  g
}
