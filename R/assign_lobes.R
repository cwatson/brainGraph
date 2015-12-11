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
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

assign_lobes <- function(g, rand=FALSE) {
  lobe <- hemi <- name <- class <- NULL
  if (!'atlas' %in% graph_attr_names(g)) {
    stop(sprintf('Input graph %s does not have an "atlas" attribute!',
                 deparse(substitute(g))))
  }
  atlas.dt <- eval(parse(text=g$atlas))

  # Check that vertex names match the atlas names
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
    if (g$atlas %in% c('dkt', 'dk', 'destrieux')) {
      counts <- atlas.dt[order(lobe), .N, by=list(lobe, hemi)]$N
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

    } else if (g$atlas %in% c('aal90', 'lpba40', 'hoa112', 'brainsuite',
                            'dk.scgm', 'dkt.scgm')) {
      counts <- atlas.dt[order(lobe), .N, by=list(lobe, hemi)]$N
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
