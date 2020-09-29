#' Count number of edges of a brain graph
#'
#' \code{count_homologous} counts the number of edges between homologous regions
#' in a brain graph (e.g. between L and R superior frontal).
#'
#' @param g A \code{brainGraph} graph object
#' @export
#' @return \code{count_homologous} - a named vector of the edge ID's connecting
#'   homologous regions
#'
#' @name Count Edges
#' @rdname count_edges
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}

count_homologous <- function(g) {
  V1 <- V2 <- NULL
  stopifnot(is.brainGraph(g))
  if (any(grepl('(\\.|_)L[0-9]*', V(g)$name)) && !grepl('destrieux', g$atlas)) {
    lh <- '(\\.|_)L[0-9]*'
    rh <- '(\\.|_)R[0-9]*'
  } else {
    lh <- '^l'
    rh <- '^r'
  }

  dt <- as.data.table(as_edgelist(g))
  eids <- dt[, which(gsub(lh, '', V1) == gsub(rh, '', V2))]
  names(eids) <- dt[eids, V1]
  return(eids)
}

#' Count number of inter-group edges
#'
#' \code{count_inter} counts the number of edges between and within all vertices
#' in one group (e.g. \emph{lobe}, \emph{hemi}, \emph{network}, etc.).
#'
#' @param group Character string specifying which grouping to calculate edge
#'   counts for. Default: \code{'lobe'}
#' @export
#' @return \code{count_inter} - a \code{data.table} of total, intra-, and
#'   inter-group edge counts
#'
#' @rdname count_edges
#' @examples
#' \dontrun{
#' g1.lobecounts <- count_inter(g[[1]][[N]], 'lobe')
#' }

count_inter <- function(g, group=c('lobe', 'hemi', 'network', 'class',
                                   'gyrus', 'Yeo_7network', 'Yeo_17network',
                                   'area', 'Brodmann')) {
  total <- intra <- inter <- NULL
  group <- match.arg(group)
  stopifnot(is.brainGraph(g), group %in% vertex_attr_names(g))

  group.names <- get(g$atlas)[, levels(get(group))]
  A <- as_adj(g, names=FALSE, sparse=FALSE)
  Nm <- length(group.names)
  mat <- matrix(0, Nm, Nm)
  vattrs <- vertex_attr(g, group)
  matches <- lapply(group.names, function(x) which(vattrs == x))
  for (i in seq_len(Nm)) {
    for (j in seq.int(i, Nm)) {
      mat[i, j] <- sum(A[matches[[i]], matches[[j]]])
    }
  }
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  intra <- diag(mat) <- diag(mat) / 2
  rownames(mat) <- colnames(mat) <- group.names
  DT <- data.table(group=group.names, intra=intra, inter=rowSums(mat)-intra)
  DT[, total := intra + inter]
  setnames(DT, 'group', group)
  return(list(mat=mat, DT=DT))
}
