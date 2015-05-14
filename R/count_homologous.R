#' Count number of edges between homologous regions of a brain graph
#'
#' This function will count the number of edges between homologous regions in a
#' brain graph (e.g. between L and R superior frontal).
#'
#' @param g An igraph graph object
#' @export
#'
#' @return An integer vector of the edge ID's connecting homologous regions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

count_homologous <- function(g) {
  atlas <- g$atlas
  eids <- vector('integer')
  if (atlas %in% c('dkt', 'dk', 'destrieux', 'brainsuite')) {
    for (i in seq_len(vcount(g) / 2)) {
      region <- substr(V(g)$name[i], 2, nchar(V(g)$name[i]))
      region.l <- paste0('l', region)
      region.r <- paste0('r', region)
      region.r.num <- which(V(g)$name == region.r)

      if (length(E(g)[i %--% region.r.num]) == 1) {
        eids <- c(eids, as.numeric(E(g)[i %--% region.r.num]))
      }
    }
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
  }
  return(eids)
}
