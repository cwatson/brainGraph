#' Count number of inter-lobar connections from a given major lobe
#'
#' This function will count the number of edges between all vertices in one
#' major lobe (e.g. Frontal) and all other major lobes.
#'
#' @param g The igraph graph object
#' @param lobe A character string indicating the lobe to count from (uppercase)
#' @export
#'
#' @return A data table of total, intra-, and inter-lobar edge counts
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g1.frontal <- count_interlobar(g1[[N]], 'Frontal')
#' }

count_interlobar <- function(g, lobe) {
  stopifnot(is_igraph(g))
  stopifnot('atlas' %in% graph_attr_names(g))

  lobe.names <- eval(parse(text=g$atlas))[, levels(lobe)]
  stopifnot(lobe %in% lobe.names)

  id <- which(lobe == lobe.names)
  total <- length(E(g)[which(V(g)$lobe == id) %--% V(g)])
  intra <- length(E(g)[which(V(g)$lobe == id) %--% which(V(g)$lobe == id)])
  inter <- total - intra

  DT <- data.table(total=total, intra=intra, inter=inter)
  return(DT)
}
