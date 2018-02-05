#' Calculate the rich club of a graph
#'
#' This function calculates the \emph{rich club} of a graph, both the
#' coefficient \eqn{\phi} and the nodes that make up this subgraph.
#'
#' @param g An \code{igraph} graph object
#' @param k Integer; the minimum degree for including a vertex (default: 1)
#' @param weighted Logical indicating whether or not edge weights should be
#'   used (default: \code{FALSE})
#' @export
#'
#' @return A list with the following components:
#' \item{phi}{The rich club coefficient, \eqn{\phi}.}
#' \item{graph}{A subgraph containing only the rich club nodes.}
#' \item{Nk}{The number of vertices in the rich club graph.}
#' \item{Ek}{The number of edges in the rich club graph.}
#'
#' @family Rich-club functions
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Zhou S., Mondragon R.J. (2004) \emph{The rich-club phenomenon
#' in the internet topology}. IEEE Comm Lett, 8:180-182.
#' @references Opsahl T., Colizza V., Panzarasa P., Ramasco J.J. (2008)
#' \emph{Prominence and control: the weighted rich-club effect}. Physical Review
#' Letters, 101.16:168702.

rich_club_coeff <- function(g, k=1, weighted=FALSE) {
  stopifnot(is_igraph(g))
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  Nv <- vcount(g)
  Nk <- sum(degs > k)
  if (Nk == 0) {
    return(list(phi=NaN, graph=make_empty_graph(), Nk=0, Ek=0))
  } else {
    rich.club.nodes <- order(degs)[(Nv - Nk + 1):Nv]
    rich.club.graph <- induced.subgraph(g, rich.club.nodes)
    Ek <- ecount(rich.club.graph)

    if (isTRUE(weighted)) {
      Wr <- sum(E(rich.club.graph)$weight)
      weights <- sort(E(g)$weight, decreasing=TRUE)[1:Ek]
      phi <- Wr / sum(weights)
    } else {
      phi <- graph.density(rich.club.graph)
    }

    return(list(phi=phi, graph=rich.club.graph, Nk=Nk, Ek=Ek))
  }
}

#' Calculate the normalized rich club coefficient
#'
#' This function will (optionally) generate a number of random graphs, calculate
#' their rich club coefficients (\eqn{\phi}), and return \eqn{\phi} of the graph
#' of interest divided by the mean across random graphs, i.e. \eqn{\phi_{norm}}.
#' If random graphs have already been generated, you can supply a list as an
#' argument (since graph generation is time consuming).
#'
#' @param g An \code{igraph} graph object
#' @param N Integer; the number of random graphs to generate (default: 100)
#' @param rand A list of \code{igraph} graph objects, if random graphs have
#'   already been generated (default: NULL)
#' @param ... Other parameters (passed to \code{\link{sim.rand.graph.par}})
#' @export
#'
#' @return A data table with columns:
#'   \item{k}{Sequence of degrees}
#'   \item{rand}{Rich-club coefficients for the random graphs}
#'   \item{orig}{Rich-club coefficients for the original graph.}
#'   \item{norm}{Normalized rich-club coefficients.}
#'   \item{p}{The P-values based on the distribution of rich-club coefficients
#'     from the random graphs.}
#'   \item{p.fdr}{The FDR-adjusted P-values}
#'   \item{density}{The observed graph's density}
#'   \item{threshold,Group,name}{(if applicable)}
#'
#' @family Rich-club functions
#' @family Random graph functions
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Colizza V., Flammini A., Serrano M.A., Vespignani A. (2006)
#' \emph{Detecting rich-club ordering in complex networks}. Nature Physics,
#' 2:110-115.

rich_club_norm <- function(g, N=1e2, rand=NULL, ...) {
  k <- orig <- p <- p.fdr <- NULL
  stopifnot(is_igraph(g))
  if (is.null(rand)) {
    rand <- sim.rand.graph.par(g, N, ...)
  } else {
    if (!all(vapply(rand, is_igraph, logical(1)))) {
      stop('Argument "rand" must be a list of igraph graph objects!')
    }
    N <- length(rand)
  }

  phi.rand <- t(sapply(rand, function(x) x$rich$phi))
  max.deg <- max(V(g)$degree)
  DT <- data.table(k=rep(seq_len(max.deg), each=N),
                   rand=c(phi.rand),
                   orig=rep(g$rich$phi, each=N))
  DT[, norm := unique(orig) / mean(rand), by=k]
  DT[, p := (sum(rand >= unique(orig)) + 1) / (N + 1), by=k]
  dt.phi <- DT[, .SD[1], by=k]
  dt.phi[, p.fdr := p.adjust(p, 'fdr')]
  dt.phi[, rand := NULL]
  dt.phi$rand <- DT[, mean(rand), by=k]$V1
  dt.phi$density <- g$density
  if ('threshold' %in% graph_attr_names(g)) dt.phi$threshold <- g$threshold
  if ('Group' %in% graph_attr_names(g)) dt.phi$Group <- g$Group
  if (!is.null(g$name)) dt.phi$Study.ID <- g$name
  return(dt.phi)
}

#' Assign graph attributes based on rich-club analysis
#'
#' This function will assign vertex- and edge-level attributes based on the
#' results of a \emph{rich-club} analysis, based on a range of vertex degrees in
#' which the rich-club coefficient was determined to be significantly greater
#' than that of a set of random graphs (see \code{\link{rich_club_norm}}).
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
#' @param deg.range Integer vector of the range of degrees indicating
#'   inclusion in the rich-club; if the default \emph{NULL}, it will be from 1
#'   to the maximum degree in the graph
#' @param adj.vsize Logical indicating whether to adjust vertex size
#'   proportional to degree (default: FALSE)
#' @export
#'
#' @return An \code{igraph} graph object with additional attributes:
#'   \item{rich}{Binary indicating membership in the rich-club}
#'   \item{type.rich}{Edge attribute indicating the type of connection}
#'   \item{color.rich}{Edge and vertex attributes}
#'   \item{size.rich}{Edge and vertex attributes}
#'
#' @family Rich-club functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' g <- rich_club_attrs(g, rich.dt[density == densities[N] & p.fdr < .01,
#'                                 range(k)])
#' }

rich_club_attrs <- function(g, deg.range=NULL, adj.vsize=FALSE) {
  stopifnot(is_igraph(g))
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

#' Calculate the rich core of a graph
#'
#' This function finds the boundary of the rich core of a graph, based on the
#' decreasing order of vertex degree. It also calculates the degree that
#' corresponds to that rank, and the core size relative to the total number of
#' vertices in the graph.
#'
#' @inheritParams rich_club_coeff
#' @export
#'
#' @return A data frame with the following components:
#' \item{density}{The density of the graph.}
#' \item{rank}{The rank of the boundary for the rich core.}
#' \item{k.r}{The degree of the vertex at the boundary.}
#' \item{core.size}{The size of the core relative to the graph size.}
#'
#' @family Rich-club functions
#'
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @references Ma A \& Mondragon R.J. (2015) \emph{Rich-cores in networks}. PLoS
#' One, 10(3): e0119678. doi:10.1371/journal.pone.0119678

rich_core <- function(g) {
  stopifnot(is_igraph(g))
  if ('degree' %in% vertex_attr_names(g)) {
    degs <- V(g)$degree
  } else {
    degs <- degree(g)
  }
  dens <- ifelse('density' %in% graph_attr_names(g), g$density, graph.density(g))
  Nv <- vcount(g)

  vorder <- order(degs, decreasing=TRUE)
  kplus <- sapply(seq_len(Nv), function(x)
                  length(E(g)[vorder[x] %--% which(degs > degs[vorder[x]])]))

  r <- max(which(kplus == max(kplus)))
  k.r <- degs[vorder][r]
  core.size <- r / Nv
  return(data.frame(density=dens, rank=r, k.r=k.r, core.size=core.size))
}
