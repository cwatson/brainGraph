#' Assign graph attributes based on rich-club analysis
#'
#' This function will assign vertex- and edge-level attributes based on the
#' results of a \emph{rich-club} analysis, based on a range of vertex degrees in
#' which the rich-club coefficient was determined to be significantly greater
#' than that of a set of random graphs (see \code{\link{rich.club.norm}}).
#'
#' Vertices which are in the rich club will be assigned an attribute
#' \code{rich}, taking on a binary value. Their colors (attribute
#' \code{color.rich}) will be either \emph{red} or \emph{gray}. Their sizes
#' (attribute \code{size.rich}) will either be 10 or will be proportional to
#' their degree.
#'
#' Edge attribute \code{type.rich} takes on three values: \emph{rich-club} (if
#' it connects two rich-club vertices), \emph{feeder} (if it connects a rich- to
#' a non-rich-club vertex), and \emph{local} (if it connects two non-rich-club
#' vertices). They will also be given a \code{color.rich} attribute (either
#' \emph{red}, \emph{orange}, or \emph{green}). Edge sizes (\code{size.rich})
#' will be largest for \emph{rich-club} connections, then smaller for
#' \emph{feeder}, and smallest for \emph{local}.
#'
#' @param g An \code{igraph} graph object
#' @param deg.range An integer vector of the range of degrees indicating
#'   inclusion in the rich-club; if the default \emph{NULL}, it will be from 1
#'   to the maximum degree in the graph
#' @param adj.vsize A logical indicating whether to adjust vertex size
#'   proportional to degree (default: FALSE)
#' @export
#'
#' @return An \code{igraph} graph object with additional attributes:
#'   \item{rich}{Binary indicating membership in the rich-club}
#'   \item{type.rich}{Edge attribute indicating the type of connection}
#'   \item{color.rich}{Edge and vertex attributes}
#'   \item{size.rich}{Edge and vertex attributes}
#' @seealso \code{\link{rich.club.norm}, \link{rich.club.coeff}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g <- rich.club.attrs(g, rich.dt[density == densities[N] & p.fdr < .01,
#'                                 range(k)])
#' }

rich.club.attrs <- function(g, deg.range=NULL, adj.vsize=FALSE) {
  if (is.null(deg.range)) deg.range <- c(1, max(V(g)$degree))
  V(g)$rich <- 0
  V(g)[V(g)$degree >= deg.range[1] & V(g)$degree <= deg.range[2]]$rich <- 1
  E(g)[which(V(g)$rich == 1) %--% which(V(g)$rich == 1)]$type.rich <- 'rich-club'
  E(g)[which(V(g)$rich == 1) %--% which(V(g)$rich == 0)]$type.rich <- 'feeder'
  E(g)[which(V(g)$rich == 0) %--% which(V(g)$rich == 0)]$type.rich <- 'local'

  V(g)[V(g)$rich == 1]$color.rich <- 'red'
  V(g)[V(g)$rich == 0]$color.rich <- 'gray'
  E(g)[E(g)$type.rich == 'rich-club']$color.rich <- 'red'
  E(g)[E(g)$type.rich == 'feeder']$color.rich <- 'orange'
  E(g)[E(g)$type.rich == 'local']$color.rich <- 'green'

  if (isTRUE(adj.vsize)) {
    V(g)[V(g)$rich == 1]$size.rich <- vec.transform(V(g)[which(V(g)$rich == 1)]$degree, 3, 15)
  } else {
    V(g)[V(g)$rich == 1]$size.rich <- 10
  }
  V(g)[V(g)$rich == 0]$size.rich <- 0
  E(g)[E(g)$type.rich == 'rich-club']$size.rich <- 3.5
  E(g)[E(g)$type.rich == 'feeder']$size.rich <- 1.5
  E(g)[E(g)$type.rich == 'local']$size.rich <- 0.5

  return(g)
}
