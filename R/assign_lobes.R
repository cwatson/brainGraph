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
#' for a random graph or a "graph of interest" (default: \code{FALSE}).
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
  stopifnot(is_igraph(g), 'atlas' %in% graph_attr_names(g))
  lobe <- hemi <- name <- index <- N <- class <- network <- x <- y <- z <- x.mni <- y.mni <- z.mni <- NULL

  DT <- eval(parse(text=g$atlas))
  nonmatches <- !V(g)$name %in% DT[, name]
  if (any(nonmatches)) {
    stop(paste('Check the following vertex names: ',
               paste(V(g)$name[nonmatches], collapse=' ')))
  }

  vorder <- match(V(g)$name, DT$name)
  lobe.nums <- DT[vorder, as.numeric(lobe)]
  V(g)$lobe <- DT[vorder, as.character(lobe)]
  V(g)$lobe.hemi <- as.numeric(DT[vorder, interaction(lobe, hemi)])
  V(g)$hemi <- DT[vorder, as.character(hemi)]

  if (isTRUE(grepl('destr', g$atlas))) V(g)$class <- DT[vorder, as.numeric(class)]
  if (g$atlas == 'dosenbach160') V(g)$network <- DT[vorder, as.character(network)]

  if (!isTRUE(rand)) {
    l.cir <- vector('integer')
    lobes <- DT[, levels(lobe)]
    V(g)$x <- V(g)$x.mni <- DT[vorder, x.mni]
    V(g)$y <- V(g)$y.mni <- DT[vorder, y.mni]
    V(g)$z <- V(g)$z.mni <- DT[vorder, z.mni]
    V(g)$color.lobe <- group.cols[lobe.nums]
    E(g)$color.lobe <- set_edge_color(g, lobe.nums)
    if (g$atlas %in% c('destrieux', 'destrieux.scgm')) {
      V(g)$color.class <- group.cols[V(g)$class]
      E(g)$color.class <- set_edge_color(g, V(g)$class)
    }
    if (g$atlas == 'dosenbach160') {
      V(g)$color.network <- group.cols[DT[vorder, as.numeric(network)]]
      E(g)$color.network <- set_edge_color(g, DT[vorder, as.numeric(network)])
      l.cir <- c(l.cir, which(V(g)$hemi == 'B'))
    }

    l.cir <- c(l.cir,
      DT[lobe == 'Frontal' & hemi == 'L', .SD[order(-y.mni, x.mni), index]],
      DT[lobe %in% c('Insula', 'Central') & hemi == 'L', .SD[order(-y.mni, x.mni), index]],
      DT[lobe %in% c('Limbic', 'Cingulate') & hemi == 'L', .SD[order(-y.mni, x.mni), index]])
    if ('SCGM' %in% lobes) {
      l.cir <- c(l.cir, DT[lobe == 'SCGM' & hemi == 'L', .SD[order(-y.mni, x.mni), index]])
    }
    l.cir <- c(l.cir,
      DT[lobe == 'Temporal' & hemi == 'L', .SD[order(-y.mni, x.mni), index]],
      DT[lobe == 'Parietal' & hemi == 'L', .SD[order(-y.mni, x.mni), index]],
      DT[lobe == 'Occipital' & hemi == 'L', .SD[order(-y.mni, x.mni), index]],
      DT[lobe == 'Occipital' & hemi == 'R', .SD[order(y.mni, x.mni), index]],
      DT[lobe == 'Parietal' & hemi == 'R', .SD[order(y.mni, x.mni), index]],
      DT[lobe == 'Temporal' & hemi == 'R', .SD[order(y.mni, x.mni), index]])
    if ('SCGM' %in% lobes) {
      l.cir <- c(l.cir, DT[lobe == 'SCGM' & hemi == 'R', .SD[order(y.mni, x.mni), index]])
    }
    l.cir <- c(l.cir,
      DT[lobe %in% c('Limbic', 'Cingulate') & hemi == 'R', .SD[order(y.mni, x.mni), index]],
      DT[lobe %in% c('Insula', 'Central') & hemi == 'R', .SD[order(y.mni, x.mni), index]],
      DT[lobe == 'Frontal' & hemi == 'R', .SD[order(y.mni, x.mni), index]])
    if ('Cerebellum' %in% lobes) {
      counts <- DT[order(lobe, hemi), .N, by=list(lobe, hemi)]
      mid1 <- counts[!lobe %in% c('Cerebellum', 'Brainstem') & hemi != 'R', sum(N)]
      mid2 <- counts[!lobe %in% c('Cerebellum', 'Brainstem') & hemi == 'R', sum(N)]
      l.cir <- c(l.cir[1:mid1],
                       which(V(g)$lobe == 'Cerebellum'),
                       l.cir[(mid1+1):(mid2+mid1)])
    }
    if ('Brainstem' %in% lobes) {
      mid1 <- counts[lobe != 'Brainstem' & hemi != 'R', sum(N)]
      mid2 <- counts[lobe != 'Brainstem' & hemi == 'R', sum(N)]
      l.cir <- c(l.cir[1:mid1],
                       which(V(g)$lobe == 'Brainstem'),
                       l.cir[(mid1+1):(mid2+mid1)])
    }
    V(g)$circle.layout <- l.cir
  }

  g
}
