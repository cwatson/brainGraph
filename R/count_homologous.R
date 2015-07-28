#' Count number of edges between homologous regions of a brain graph
#'
#' This function will count the number of edges between homologous regions in a
#' brain graph (e.g. between L and R superior frontal).
#'
#' @param g An igraph graph object
#' @export
#'
#' @return A named vector of the edge ID's connecting homologous regions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

count_homologous <- function(g) {
  if (!'atlas' %in% graph_attr_names(g)) {
    stop(sprintf('Graph "%s" must have an "atlas" attribute!',
                 deparse(substitute(g))))
  }
  atlas <- g$atlas
  Nv <- vcount(g)
  if (atlas %in% c('dkt', 'dk', 'destrieux', 'brainsuite')) {
    verts <- seq_len(Nv / 2)
    regions <- substr(V(g)[verts]$name, 2, nchar(V(g)[verts]$name))
    regions.l <- paste0('l', regions)
    regions.r <- paste0('r', regions)

    eids <- unlist(Map(function(x, y)
                          as.numeric(E(g)[x %--% y]),
                          regions.l,
                          regions.r))

  } else if (atlas %in% c('aal90', 'hoa112')) {
    for (i in seq(1, vcount(g), by=2)) {
      if (length(E(g)[i %--% (i+1)]) == 1) {
        eids <- c(eids, as.numeric(E(g)[i %--% (i+1)]))
      }
    }
  } else if (atlas == 'lpba40') {
    for (i in seq(1, vcount(g) - 2, by=2)) {
      if (length(E(g)[i %--% (i+1)]) == 1) {
        eids <- c(eids, as.numeric(E(g)[i %--% (i+1)]))
      }
    }
  } else {
    stop(sprintf('Atlas "%s" is not a valid name!', atlas))
  }
  return(eids)
}
