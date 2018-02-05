#' Save PNG of three views of a brain graph
#'
#' This function will save a PNG file to disk containing three views (columns)
#' of a brain graph (from left-to-right): left sagittal, axial, and right
#' sagittal. The number of rows in the figure will equal the number of groups to
#' plot.
#'
#' The function argument \code{N} tells the function to use the \code{N-th}
#' element of the input list \code{g.list} for each group. So, for example, if
#' \code{g.list} consists of lists of graphs for two groups, and \code{N} is 4,
#' then the plots for \code{g.list[[1]][[4]]} and \code{g.list[[2]][[4]]} will
#' be written to the file.
#'
#' The \code{subgraph} argument can be used to apply one or more conditions for
#' subsetting the graph. If you would like multiple conditions, then it must be
#' a \code{list} variable that equals in length to the number of groups. For a
#' single group and multiple conditions, simply write e.g., \code{groups=c(1,
#' 1)}. The \code{main} argument has the same rule except it controls the main
#' plot title, which appears in the \emph{axial} view along with the
#' \emph{Group} name.
#'
#' @param g.list A list of lists of \code{igraph} graph objects
#' @param groups An integer vector indicating which groups to plot; corresponds
#'   to the first element of the list \code{g.list} (default: 1)
#' @param N Integer corresponding to the second element of the list
#'   \code{g.list} (default: 1)
#' @param filename Character string of the filename of the PNG to be written
#'   (default: 'tmp.png')
#' @param subgraph A list of character strings to (optionally) subset the
#'   graph(s), possibly by multiple conditions (default: \code{NULL})
#' @param main A list of character strings to be placed in the main title of the
#'   center plot for each group (default: \code{NULL})
#' @param ... Other arguments passed to
#'   \code{\link{plot_brainGraph}}
#' @export
#' @importFrom oro.nifti nifti
#'
#' @family Plotting functions
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' plot_brainGraph_multi(g.hubs, groups=1:2, filename='Figure01_hubs.png',
#'   subgraph='N > 0', vertex.color='color.lobe', vertex.size=15,
#'   show.legend=TRUE, vertex.label.cex=1.5)
#' ## Single group, different subgraphs for each plot
#' plot_brainGraph_multi(g, groups=c(1, 1), N=5, filename='5_6core.png',
#'   vertex.color='color.lobe', edge.color='color.lobe', vertex.label=NA,
#'   subgraph=list('coreness > 5', 'coreness > 6'),
#'   main=list('k-core 5', 'k-core 6'))
#' }

plot_brainGraph_multi <- function(g.list, groups=1, N=1, filename='tmp.png',
                                  subgraph=NULL, main=NULL, ...) {
  stopifnot('Group' %in% graph_attr_names(g.list[[groups[1]]][[N]]))

  X <- mni152@.Data
  L <- nifti(X[rev(seq_len(nrow(X))), rev(seq_len(ncol(X))), ])
  R <- aperm(X, c(2, 3, 1))
  L <- aperm(L, c(2, 3, 1))

  zlim <- c(3500, max(X[, , 46]))
  imcol <- gray(0:64/64)
  breaks <- c(zlim[1], seq(min(zlim), max(zlim), length=64), zlim[2])

  kNumGroups <- length(groups)
  png(filename=filename, width=24, height=8*kNumGroups, units='in', res=300)
  layout(matrix(1:(3*kNumGroups), kNumGroups, 3, byrow=TRUE),
         widths=c(4/3, 1, 4/3), heights=lcm(8*2.54))

  # See if there are multiple "subgraph" or "main" arguments
  if (!is.null(subgraph)) {
    if (kNumGroups > 1) {
      if (length(subgraph) > 1) {
        stopifnot(is.list(subgraph),
                  length(subgraph) == kNumGroups)
      } else if (length(subgraph) == 1) {
        subgraph <- rep(as.list(subgraph), kNumGroups)
      }
    } else {
      subgraph <- as.list(subgraph)
    }
  } else {
    subgraph <- vector('list', length=kNumGroups)
  }
  if (!is.null(main)) {
    if (kNumGroups > 1) {
      if (length(main) > 1) {
        stopifnot(is.list(main),
                  length(main) == kNumGroups)
      } else if (length(main) == 1) {
        main <- rep(as.list(main), kNumGroups)
      }
    } else {
      main <- as.list(main)
    }
  } else {
    main <- NULL
  }

  for (i in seq_along(groups)) {
    main.title <- paste0('\n', g.list[[groups[i]]][[N]]$Group)
    if (!is.null(main)) {
      main.title <- paste0(main.title, ': ', main[[i]])
    }
    # Left sag.
    par(mar=c(0, 0, 0, 0), bg='black')
    graphics::image(1:109, 1:91, L[, , 30], col=imcol, breaks=breaks, asp=0)
    par(new=TRUE)
    plot(g.list[[groups[i]]][[N]], plane='sagittal', hemi='L', main='\n\n\nLH',
                    subgraph=subgraph[[i]], ...)

    # Axial
    par(mar=c(0, 0, 0, 0), bg='black')
    graphics::image(1:91, 1:109, X[, , 46], col=imcol, breaks=breaks)
    par(new=TRUE)
    plot(g.list[[groups[i]]][[N]], main=main.title, subgraph=subgraph[[i]], ...)

    # Right sag.
    par(mar=c(0, 0, 0, 0), bg='black')
    graphics::image(1:109, 1:91, R[, , 30], col=imcol, breaks=breaks)
    par(new=TRUE)
    plot(g.list[[groups[i]]][[N]], plane='sagittal', hemi='R', main='\n\n\nRH',
                    subgraph=subgraph[[i]], ...)
  }
  dev.off()
}
