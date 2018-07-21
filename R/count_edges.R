#' Count number of edges of a brain graph
#'
#' \code{count_homologous} counts the number of edges between homologous regions
#' in a brain graph (e.g. between L and R superior frontal).
#'
#' @param g An \code{igraph} graph object
#' @export
#'
#' @return \code{count_homologous} - a named vector of the edge ID's connecting
#'   homologous regions
#'
#' @name CountEdges
#' @aliases count_homologous
#' @rdname count_edges
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

count_homologous <- function(g) {
  stopifnot(inherits(g, 'brainGraph'))

  eids <- unlist(Map(function(x, y)
                     as.numeric(E(g)[x %--% y]),
                     which(V(g)$hemi == 'L'),
                     which(V(g)$hemi == 'R')))
  names(eids) <- as_edgelist(g)[eids]
  return(eids)
}

#' Count number of inter-lobar connections from a given major lobe
#'
#' \code{count_interlobar} counts the number of edges between all vertices in
#' one major lobe (e.g. Frontal) and all other major lobes.
#'
#' @param lobe A character string indicating the lobe to count from (uppercase)
#' @export
#'
#' @return \code{count_interlobar} - a \code{data.table} of total, intra-, and
#'   inter-lobar edge counts
#'
#' @aliases count_interlobar
#' @rdname count_edges
#' @examples
#' \dontrun{
#' g1.frontal <- count_interlobar(g[[1]][[N]], 'Frontal')
#' }

count_interlobar <- function(g, lobe) {
  stopifnot(inherits(g, 'brainGraph'))

  lobe.names <- get(g$atlas)[, levels(lobe)]
  stopifnot(lobe %in% lobe.names)

  total <- length(E(g)[which(V(g)$lobe == lobe) %--% V(g)])
  intra <- length(E(g)[which(V(g)$lobe == lobe) %--% which(V(g)$lobe == lobe)])
  inter <- total - intra

  DT <- data.table(total=total, intra=intra, inter=inter)
  return(DT)
}
