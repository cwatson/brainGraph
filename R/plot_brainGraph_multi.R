#' Save PNG of one or three views for all graphs in a brainGraphList
#'
#' \code{plot_brainGraph_multi} writes a PNG file to disk containing three views
#' (columns) of 1 or more \code{brainGraph} objects (from left-to-right): left
#' sagittal, axial, and right sagittal. The number of rows in the figure will
#' equal the number of graphs to plot.
#'
#' Whether the first input is a \code{brainGraphList} object or a \code{list} of
#' \code{brainGraph} objects, \emph{all} graphs in the object will be displayed
#' in the figure. For \code{plot_brainGraph_multi}, this may be undesirable if
#' you have more than 4 or 5 graphs in one object. You can choose fewer by using
#' simple subsetting operations (see Examples below).
#'
#' @section Using subgraphs, titles, and labels:
#' There are three arguments that can differ for each graph to be displayed.
#' Each follows the same \dQuote{rules}. If you would like the same value
#' applied to all graphs, you can specify a \code{character} string. If you
#' would like a different value for each group, you must supply a \code{vector}
#' or \code{list} with length equal to the number of graphs. If its length is
#' less than the number of graphs, values will be recycled. To \dQuote{skip}
#' applying a value to one (or more) graph(s), you can use the \code{NULL} value
#' only within a list (see the Examples below).
#' \describe{
#'   \item{subgraph}{Can be used to apply one or more conditions for subsetting
#'     the graph(s).}
#'   \item{main}{Controls the main plot title, which appears in the \emph{axial}
#'     view along with each graph's \code{name} attribute. Depending on the
#'     \code{level} of the \code{brainGraphList}, this will either be a Study
#'     ID, Group name, or contrast name.}
#'   \item{label}{Can be used to print a text label in a corner for each
#'     group/graph. For example, you can print a letter if you will refer to,
#'     e.g., \dQuote{Figure 1A}, \dQuote{Figure 1B}, etc.}
#' }
#'
#' @note All other arguments (passed to \code{\link{plot.brainGraph}}) will be
#' applied to \emph{all} graphs. For example, if you include
#' \code{vertex.label=NA} in the function call, vertex labels will be omitted
#' for all graphs.
#'
#' @param g.list A \code{brainGraphList} or a list of \code{brainGraph} objects
#' @param filename Character string of the filename of the PNG to be written.
#' @param subgraph A vector or list of character strings to (optionally) subset
#'   the graph(s), possibly by multiple conditions
#' @param main A vector or list of character strings to be placed in the main
#'   title of the center (axial) plot for each graph
#' @param label A vector or list of character strings to be placed in one of the
#'   corners of the left plot (sagittal) in each row
#' @param cex.main Numeric specifying the level of character expansion for the
#'   plot titles. Default: \code{1} (no expansion)
#' @param ... Other arguments passed to \code{\link{plot.brainGraph}}
#' @export
#' @importFrom graphics layout par
#' @importFrom grDevices dev.off png
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' ## "g.hubs" contains 2 groups; apply same subset to both
#' plot_brainGraph_multi(g.hubs, filename='Figure01_hubs.png',
#'   subgraph='N > 0', vertex.color='color.lobe', vertex.size=15,
#'   show.legend=TRUE, vertex.label.cex=1.5)
#'
#' ## Single group, different subgraphs for both plots
#' ## "g" is a "brainGraphList" object
#' gg <- g[rep(1, 3), drop=FALSE]
#' plot_brainGraph_multi(gg, filename='group1_5-6-7core.png',
#'   vertex.color='color.lobe', edge.color='color.lobe', vertex.label=NA,
#'   subgraph=as.list(paste('coreness >', 5:7)),
#'   main=as.list(paste('k-core', 5:7)))
#'
#' ## Apply different subset for groups 1 & 3; no subset for group 2
#' plot_brainGraph_multi(g, groups=1:3, vertex.label=NA,
#'   subgraph=list('degree > 5', NULL, 'degree > 4'))
#' }

plot_brainGraph_multi <- function(g.list, filename='orthoview.png',
                                  subgraph=NULL, main=NULL, label=NULL,
                                  cex.main=1, ...) {
  if (!inherits(g.list, 'brainGraphList')) try(g.list <- as_brainGraphList(g.list))
  g.list <- g.list[]
  kNumGraphs <- length(g.list)

  width <- 671 * 3
  height <- 673
  png(filename=filename, width=width, height=height*kNumGraphs, units='px', res=300)
  layout(matrix(1:(3*kNumGraphs), kNumGraphs, 3, byrow=TRUE),
         widths=c(4/3, 1, 4/3), heights=rep.int(1, kNumGraphs))

  # See if there are multiple "subgraph", "main", or "label" arguments
  if (!is.null(subgraph)) subgraph <- rep_len(as.list(subgraph), kNumGraphs)
  main.title <- lapply(g.list, graph_attr, 'name')
  if (!is.null(main)) {
    main <- rep_len(as.list(main), kNumGraphs)
    main.title <- Map(paste0, main.title, ': ', main)
  }
  if (!is.null(label)) label <- rep_len(as.list(label), kNumGraphs)

  for (i in seq_len(kNumGraphs)) {
    plot(g.list[[i]], plane='sagittal', hemi='L', subgraph=subgraph[[i]], main='\n\n\nLH',
         label=label[[i]], cex.main=cex.main, ...)
    plot(g.list[[i]], subgraph=subgraph[[i]], main=main.title[[i]], cex.main=cex.main, ...)
    plot(g.list[[i]], plane='sagittal', hemi='R', subgraph=subgraph[[i]], main='\n\n\nRH',
         cex.main=cex.main, ...)
  }
  dev.off()
}

#' Save PNG of a single view for all graphs in a brainGraphList
#'
#' \code{slicer} writes a PNG file to disk containing a single view (i.e.,
#' either sagittal, axial, or circular) of all \code{brainGraph} objects in the
#' input list/\code{brainGraphList}.
#'
#' @param nrows Integer; the number of rows in the figure
#' @param ncols Integer; the number of columns in the figure
#' @inheritParams plot.brainGraph
#' @export
#' @rdname plot_brainGraph_multi

slicer <- function(g.list, plane, hemi, nrows, ncols, filename='all.png',
                   main=NULL, cex.main=1, ...) {
  if (!inherits(g.list, 'brainGraphList')) try(g.list <- as_brainGraphList(g.list))
  g.list <- g.list[]
  kNumGraphs <- length(g.list)
  kNumPlots <- nrows * ncols
  if (kNumPlots < kNumGraphs) {
    warning('Specified # of rows and columns less than the # of graphs; adjusting nrows.')
    rem <- (kNumGraphs / nrows) %% ncols
    nrows <- (kNumGraphs %/% ncols) + (rem > 0)
    kNumPlots <- nrows * ncols
  }

  width <- 671
  height <- if (plane == 'axial') 768 else 640
  png(filename=filename, width=width*ncols, height=height*nrows, units='px', res=300)
  layout(matrix(seq_len(kNumPlots), nrows, ncols, byrow=TRUE),
         widths=rep.int(1, ncols), heights=rep.int(1, nrows))

  main.title <- lapply(g.list, graph_attr, 'name')
  if (!is.null(main)) {
    main <- rep_len(as.list(main), kNumGraphs)
    main.title <- Map(paste0, main.title, ': ', main)
  }
  for (i in seq_len(kNumGraphs)) {
    plot(g.list[[i]], plane=plane, hemi=hemi, main=main.title[[i]], cex.main=cex.main, ...)
  }
  if (kNumPlots > kNumGraphs) {
    for (i in seq_len(kNumPlots - kNumGraphs)) par(new=TRUE)
  }
  dev.off()
}
