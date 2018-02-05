#' Create a data table with graph global and vertex measures
#'
#' \code{graph_attr_dt} is a helper function that takes a list of graphs and
#' creates a \code{data.table} of global measures for each graph. Each row will
#' be for a different graph.
#'
#' @param g.list A list of \code{igraph} graph objects
#' @param group A character string indicating group membership (default:
#'   \code{NULL})
#' @export
#'
#' @return A \code{data.table}
#'
#' @name DataTables
#' @aliases graph_attr_dt
#' @rdname data_tables
#'
#' @seealso \code{\link[igraph]{graph_attr}, \link[igraph]{graph_attr_names}}

graph_attr_dt <- function(g.list, group=NULL) {
  Group <- NULL
  inds <- which(sapply(graph_attr(g.list[[1]]), class) %in% c('numeric', 'integer'))
  g.attrs <- graph_attr_names(g.list[[1]])
  g.attr.num <- g.attrs[inds]
  g.dt <- as.data.table(sapply(g.attr.num, function(x)
                               sapply(g.list, graph_attr, x)))

  if (length(g.list) == 1) {
    g.dt <- as.data.table(t(g.dt))
    colnames(g.dt) <- g.attr.num
  }

  if ('name' %in% g.attrs) g.dt$Study.ID <- sapply(g.list, function(x) x$name)
  if ('atlas' %in% g.attrs) g.dt$atlas <- g.list[[1]]$atlas
  if ('modality' %in% g.attrs) g.dt$modality <- sapply(g.list, function(x) x$modality)
  if ('weighting' %in% g.attrs) g.dt$weighting <- g.list[[1]]$weighting
  if (is.null(group)) {
    if ('Group' %in% g.attrs) g.dt$Group <- sapply(g.list, function(x) x$Group)
  } else {
    g.dt$Group <- group
  }
  return(g.dt)
}

#' Create a data table with graph vertex measures
#'
#' \code{vertex_attr_dt} is a helper function that creates a \code{data.table}
#' in which each row is a vertex and each column is a different network measure
#' (degree, centrality, etc.). It is partly a wrapper for
#' \code{\link[igraph]{as_data_frame}}.
#'
#' @param g An \code{igraph} graph object
#' @export
#'
#' @aliases vertex_attr_dt
#' @rdname data_tables
#'
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names},
#' \link[igraph]{as_data_frame}}
#' @examples
#' \dontrun{
#' dt.V <- vertex_attr_dt(g)
#' setcolorder(dt.V, c('modality', 'atlas', 'Group', names(dt.V)[1:28]))
#' }

vertex_attr_dt <- function(g, group=NULL) {
  lobe <- name <- Group <- network <- NULL
  atlas.dt <- get(g$atlas)

  dt.V <- setDT(as_data_frame(g, what='vertices'))
  cols.char <- names(which(sapply(vertex_attr(g), class) == 'character'))
  cols.rem <- c(cols.char, 'x', 'y', 'z', 'x.mni', 'y.mni', 'z.mni',
                'lobe.hemi', 'lobe', 'circle.layout', 'comm', 'comp', 'circle.layout.comm')
  if (is_weighted(g)) cols.rem <- c(cols.rem, 'comm.wt')
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
