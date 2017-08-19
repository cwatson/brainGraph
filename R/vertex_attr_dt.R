#' Create a data table with graph vertex measures
#'
#' This is just a helper function that creates a \code{data.table} in which each
#' row is a vertex and each column is a different network measure (degree,
#' centrality, etc.). It is partly a wrapper for
#' \code{\link[igraph]{as_data_frame}}.
#'
#' @param g An \code{igraph} graph object
#' @param group A character string indicating group membership (default:
#'   \code{NULL})
#' @export
#'
#' @return A data table; each row is for a different vertex
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names},
#' \link[igraph]{as_data_frame}}
#' @examples
#' \dontrun{
#' dt.V <- vertex_attr_dt(g)
#' setcolorder(dt.V, c('modality', 'atlas', 'Group', names(dt.V)[1:28]))
#' }

vertex_attr_dt <- function(g, group=NULL) {
  lobe <- name <- Group <- network <- NULL
  atlas.dt <- eval(parse(text=g$atlas))

  dt.V <- setDT(as_data_frame(g, what='vertices'))
  cols.char <- names(which(sapply(vertex_attr(g), class) == 'character'))
  cols.rem <- c(cols.char, 'x', 'y', 'z', 'x.mni', 'y.mni', 'z.mni',
                'lobe.hemi', 'lobe', 'circle.layout', 'comm', 'comm.wt', 'comp',
                'circle.layout.comm')
  cols.rem <- setdiff(cols.rem, c('name', 'hemi'))
  dt.V[, eval(cols.rem) := NULL]
  if (isTRUE(grepl('destr', g$atlas))) {
    dt.V$class <- atlas.dt[, levels(class)][V(g)$class]
  }
  if ('network' %in% cols.char) dt.V$network <- V(g)$network
  dt.V$density <- g$density
  dt.V$lobe <- V(g)$lobe
  setnames(dt.V, 'name', 'region')
  setcolorder(dt.V,
              c('density', 'region', 'lobe', 'hemi',
                names(dt.V[, !c('density', 'region', 'lobe', 'hemi'), with=F])))

  if ('name' %in% graph_attr_names(g)) dt.V$Study.ID <- g$name
  if ('modality' %in% graph_attr_names(g)) dt.V$modality <- g$modality
  if ('weighting' %in% graph_attr_names(g)) dt.V$weighting <- g$weighting
  if ('threshold' %in% graph_attr_names(g)) dt.V$threshold <- g$threshold
  if ('atlas' %in% graph_attr_names(g)) dt.V$atlas <- g$atlas
  if (is.null(group)) {
    if ('Group' %in% graph_attr_names(g)) dt.V$Group <- g$Group
  } else {
    dt.V$Group <- group
  }
  return(dt.V)
}
