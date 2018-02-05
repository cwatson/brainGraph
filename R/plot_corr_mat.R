#' Plot a correlation matrix
#'
#' This function will plot a correlation matrix in the form of a ``heatmap''.
#' You have the choice to plot the vertices in an order based on either
#' community or lobe membership, and they will be colored accordingly.
#'
#' @param corrs The correlation matrix
#' @param ordered A logical indicating whether or not to order vertices
#' (default:TRUE)
#' @param type Character string, one of: 'comm', 'comm.wt', 'lobe', or 'network'
#' @param g An igraph graph object; not required if \emph{ordered} is FALSE
#' @param group A character vector of the group name (default: NULL)
#' @export
#'
#' @return A ggplot object
#' @seealso \code{\link[ggplot2]{geom_tile}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' matplot1 <- plot_corr_mat(corrs[[1]]$r.thresh[, , N], g=g[[1]][[N]],
#'                           group=groups[1])
#' }

plot_corr_mat <- function(corrs, ordered=TRUE,
                          type=c('comm', 'comm.wt', 'lobe', 'network'),
                          g=NULL, group=NULL) {
  Var1 <- Var2 <- memb <- color <- color.test <- legend.t <- value <- color.text <- NULL
  base_size <- ifelse(nrow(corrs) > 90, 7.5, 9)

  if (isTRUE(ordered)) {
    stopifnot(!is.null(g))
    if (is.null(rownames(corrs))) rownames(corrs) <- colnames(corrs) <- V(g)$name
    Nv <- nrow(corrs)
    cols <- group.cols

    type <- match.arg(type)
    create.dt <- function(dat, graph, v.attr) {
      lobe <- memb1 <- memb2 <- Var1 <- Var2 <- value <- legend.t <- color <- color.text <- NULL
      if (v.attr %in% c('comm', 'comm.wt')) {
        memb <- vertex_attr(graph, v.attr)
        tab <- table(memb)
        group.nums <- as.integer(names(tab))
        group.max <- length(group.nums)
        group.nums <- c(group.nums, group.max + c(1, 2))
        new.order <- order(match(memb, group.nums))
        legend.title <- 'Communities (#)'
      } else if (v.attr %in% c('lobe', 'network')) {
        atlas.dt <- get(graph$atlas)
        memb <- atlas.dt[, as.numeric(get(v.attr))]
        group.nums <- c(atlas.dt[, levels(get(v.attr))])
        group.max <- length(group.nums)
        group.nums <- c(group.nums, 'Inter', '')
        new.order <- order(memb)
        legend.title <- tools::toTitleCase(v.attr)
      }
      names(memb) <- V(graph)$name
      dat <- dat[new.order, new.order]
      dat.m <- melt(dat)
      setDT(dat.m)
      dat.m[, memb1 := group.nums[memb[as.character(Var1)]]]
      dat.m[, memb2 := group.nums[memb[as.character(Var2)]]]
      dat.m[, memb := ifelse(value == 1,
                             ifelse(memb1 == memb2, memb1, group.nums[group.max + 1]),
                             group.nums[group.max + 2])]
      if (v.attr %in% c('comm', 'comm.wt')) {
        dat.m[, memb := factor(memb, levels=seq_len(max(memb)))]
      } else {
        dat.m[, memb1 := factor(memb1, levels=group.nums)]
        dat.m[, memb := factor(memb, levels=group.nums)]
      }
      dat.m[, legend.t := legend.title]
      cols.new <- c(cols[seq_len(nlevels(dat.m$memb) - 2)], 'gray50', 'white')
      dat.m[, color := cols.new[as.numeric(memb)]]
      dat.m[, color := factor(color, levels=cols.new)]
      dat.m[, color.text := cols[as.numeric(memb1)]]
      return(dat.m)
    }
    corrs.m <- create.dt(corrs, g, type)

    matplot <- ggplot(corrs.m, aes(Var1, Var2, fill=memb)) +
      geom_tile() +
      scale_fill_manual(values=corrs.m[, levels(color)]) +
      theme(axis.ticks=element_blank(),
           axis.text.x=element_text(size=0.7*base_size, angle=45,
                                    color=corrs.m[1:Nv, color.text]),
           axis.title.x=element_blank(),
           axis.text.y=element_text(size=0.7*base_size,
                                    color=corrs.m[rev(1:Nv), color.text]),
           axis.title.y=element_blank(),
           plot.title=element_text(hjust=0.5, face='bold')) +
      labs(title=group, fill=corrs.m[, unique(legend.t)]) +
      ylim(rev(levels(corrs.m$Var2)))

  } else {
    corrs.m <- melt(corrs)
    matplot <- ggplot(corrs.m, aes(Var1, Var2, fill=value)) +
      geom_tile() +
      theme(axis.ticks=element_blank(),
           axis.text.x=element_text(size=0.7*base_size, angle=45),
           axis.title.x=element_blank(),
           axis.text.y=element_text(size=0.7*base_size),
           axis.title.y=element_blank(),
           plot.title=element_text(hjust=0.5, face='bold')) +
      labs(title=group) +
      ylim(rev(levels(corrs.m$Var2)))
    if (identical(sum(abs(corrs)) - sum(corrs == 1), 0)) {
      matplot <- matplot + scale_fill_gradient2(low='white', high='blue') +
        theme(legend.position='none')
    }
  }
  return(matplot)
}
