#' Create a data table with graph global and vertex measures
#'
#' \code{graph_attr_dt} is a helper function that takes a \code{brainGraphList}
#' or a list of graphs and creates a \code{data.table} of global measures for
#' each graph. Each row will be for a different graph.
#'
#' @param bg.list A \code{brainGraphList} object, or a list of graph objects
#' @export
#' @return A \code{data.table}
#'
#' @name DataTables
#' @aliases graph_attr_dt
#' @rdname data_tables
#'
#' @seealso \code{\link[igraph]{graph_attr}, \link[igraph]{graph_attr_names}}

graph_attr_dt <- function(bg.list) {
  name <- NULL
  if (inherits(bg.list, 'brainGraphList')) {
    level <- bg.list$level
    bg.list <- bg.list$graphs
  } else {
    if (!inherits(bg.list, 'list')) bg.list <- list(bg.list)
  }
  N <- length(bg.list)
  Group <- NULL
  inds <- which(vapply(graph_attr(bg.list[[1]]), is.numeric, logical(1)))
  g.attrs <- graph_attr_names(bg.list[[1]])
  g.attr.num <- names(inds)
  g.dt <- as.data.table(vapply(g.attr.num, function(x)
                               vapply(bg.list, graph_attr, numeric(1), x), numeric(N)))

  if (N == 1) {
    g.dt <- as.data.table(t(g.dt))
    colnames(g.dt) <- g.attr.num
  }

  g.attrs.char <- c('name', 'atlas', 'modality', 'weighting', 'Group')
  for (x in g.attrs.char) {
    if (x %in% g.attrs) g.dt[, eval(x) := vapply(bg.list, graph_attr, character(1), x)]
  }
  if (level == 'subject') {
    setnames(g.dt, 'name', 'Study.ID')
  } else if (level == 'group') {
    g.dt[, name := NULL]
  }

  return(g.dt)
}

#' Create a data table with graph vertex measures
#'
#' \code{vertex_attr_dt} is a helper function that creates a \code{data.table}
#' in which each row is a vertex and each column is a different network measure
#' (degree, centrality, etc.).
#'
#' @inheritParams graph_attr_dt
#' @export
#'
#' @aliases vertex_attr_dt
#' @rdname data_tables
#'
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names},
#' \link[igraph]{graph_from_data_frame}}
#' @examples
#' \dontrun{
#' dt.V <- vertex_attr_dt(g)
#' setcolorder(dt.V, c('modality', 'atlas', 'Group', names(dt.V)[1:28]))
#' }

vertex_attr_dt <- function(bg.list) {
  lobe <- name <- Group <- network <- NULL
  if (inherits(bg.list, 'brainGraphList')) {
    level <- bg.list$level
    bg.list <- bg.list$graphs
    atlas <- bg.list$atlas
  } else {
    if (!inherits(bg.list, 'list')) bg.list <- list(bg.list)
    atlas <- bg.list[[1]]$atlas
    level <- 'subject'
  }

  dt.V <- rbindlist(lapply(bg.list, as_data_frame, what='vertices'))
  cols.char <- names(which(sapply(vertex_attr(bg.list[[1]]), is.character)))
  cols.rem <- setdiff(cols.char, c('name', 'lobe', 'hemi', 'class', 'network'))
  cols.rem <- c(cols.rem, 'x', 'y', 'z', 'x.mni', 'y.mni', 'z.mni',
                'lobe.hemi', 'circle.layout', 'circle.layout.comm')
  dt.V[, eval(cols.rem) := NULL]
  setnames(dt.V, 'name', 'region')

  # Add some important graph attributes, as well
  g.attrs <- graph_attr_names(bg.list[[1]])
  g.attrs.char <- c('name', 'atlas', 'modality', 'weighting', 'Group', 'density', 'threshold')
  for (x in g.attrs.char) {
    if (x %in% g.attrs) dt.V[, eval(x) := sapply(bg.list, graph_attr, x)]
  }
  if (level == 'subject') {
    setnames(dt.V, 'name', 'Study.ID')
  } else if (level == 'group') {
    dt.V[, name := NULL]
  }
  setcolorder(dt.V,
              c('density', 'region', 'lobe', 'hemi',
                names(dt.V[, !c('density', 'region', 'lobe', 'hemi'), with=F])))

  return(dt.V)
}
