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
  if (!is.igraph(g)) {
    stop(sprintf('%s is not a graph object', deparse(substitute(g))))
  }
  if (!'atlas' %in% graph_attr_names(g)) {
    stop(sprintf('Graph "%s" must have an "atlas" attribute',
                 deparse(substitute(g))))
  }
  atlas <- g$atlas
  Nv <- vcount(g)
  if (atlas %in% c('dkt', 'dk', 'destrieux', 'brainsuite',
                   'dk.scgm', 'dkt.scgm')) {
    half <- Nv / 2
    regions.l <- V(g)$name[seq_len(half)]
    regions.r <- V(g)$name[(half + 1):Nv]

  } else if (atlas %in% c('aal90', 'hoa112')) {
    regions.l <- V(g)$name[seq(1, Nv, by=2)]
    regions.r <- V(g)$name[seq(2, Nv, by=2)]

  } else if (atlas == 'lpba40') {
    regions.l <- V(g)$name[seq(1, Nv - 2, by=2)]
    regions.r <- V(g)$name[seq(2, Nv - 2, by=2)]

  } else {
    stop(sprintf('Atlas "%s" is not a valid name', atlas))
  }

  eids <- unlist(Map(function(x, y)
                        as.numeric(E(g)[x %--% y]),
                        regions.l,
                        regions.r))
  return(eids)
}
