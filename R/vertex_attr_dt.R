#' Create a data table with graph vertex measures
#'
#' This is just a helper function that creates a data table in which each row is
#' a vertex and each column is a different network measure (degree, centrality,
#' etc.).
#'
#' @param g An \code{igraph} graph object
#' @param group A character string indicating group membership (default: NULL)
#' @export
#'
#' @return A data table; each row is for a different vertex
#' @seealso \code{\link[igraph]{vertex_attr}, \link[igraph]{vertex_attr_names},
#' \link[igraph]{as_data_frame}}

vertex_attr_dt <- function(g, group=NULL) {
  lobe <- name <- Group <- NULL
  atlas.dt <- eval(parse(text=g$atlas))

  net.meas <- setDT(as_data_frame(g, what='vertices'))
  net.meas[, c('x', 'y', 'z', 'x.mni', 'y.mni', 'z.mni', 'lobe.hemi',
               'circle.layout', 'comm', 'comm.wt', 'comp', 'circle.layout.comm',
               'color.comm', 'color.comm.wt', 'color.comp', 'color.lobe') := NULL]
  if (g$atlas == 'destrieux') {
    net.meas[, 'color.class' := NULL]
    net.meas$class <- atlas.dt[, levels(class)][V(g)$class]
  }
  net.meas$density <- g$density
  net.meas$lobe <- atlas.dt[, levels(lobe)][V(g)$lobe]
  setnames(net.meas, 'name', 'region')
  setcolorder(net.meas,
              c('density', 'region', 'lobe', 'hemi',
                names(net.meas[, !c('density', 'region', 'lobe', 'hemi'),
                      with=F])))

  if ('name' %in% graph_attr_names(g)) net.meas$Study.ID <- g$name
  if ('modality' %in% graph_attr_names(g)) net.meas$modality <- g$modality
  if ('atlas' %in% graph_attr_names(g)) net.meas$atlas <- g$atlas
  if (is.null(group)) {
    if ('Group' %in% graph_attr_names(g)) {
      net.meas$Group <- g$Group
      setkey(net.meas, 'region', 'lobe', 'hemi', Group)
    }
    setkey(net.meas, 'region', 'lobe', 'hemi')
  } else {
    net.meas$Group <- group
    setkey(net.meas, 'region', 'lobe', 'hemi', Group)
  }
  return(net.meas)
}
