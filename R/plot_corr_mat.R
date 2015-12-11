#' Plot a correlation matrix
#'
#' This function will plot a correlation matrix in the form of a ``heatmap''.
#' You have the choice to plot the vertices in an order based on either
#' community or lobe membership, and they will be colored accordingly.
#'
#' @param corrs The correlation matrix
#' @param ordered A logical indicating whether or not to order vertices
#' (default:TRUE)
#' @param type Character string, either 'comm' or 'lobe'
#' @param g An igraph graph object; not required if \emph{ordered} is FALSE
#' @param group A character vector of the group name (default: NULL)
#' @export
#'
#' @return A ggplot object
#' @seealso \code{\link[ggplot2]{geom_tile}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' matplot1 <- plot_corr_mat(corrs[[1]][[N]]$r.thresh, g=g[[1]][[N]],
#'                            group=groups[1])
#' }

plot_corr_mat <- function(corrs, ordered=TRUE, type=c('comm', 'lobe'), g=NULL,
                          group=NULL) {
  Var1 <- Var2 <- memb <- color <- color.test <- legend.t <- value <- color.text <- NULL
  base_size <- ifelse(nrow(corrs) > 90, 7.5, 9)

  if (isTRUE(ordered)) {
    if (is.null(g)) {
      stop('You must provide a graph object for vertex ordering!')
    }
    Nv <- nrow(corrs)
    cols <- group.cols

    type <- match.arg(type)
    create.dt <- function(dat, graph, v.attr) {
      lobe <- memb1 <- memb2 <- Var1 <- Var2 <- value <- legend.t <- color <- color.text <- NULL
      memb <- vertex_attr(graph, v.attr)
      if (v.attr == 'comm') {
        tab <- table(memb)
        group.nums <- as.integer(names(tab))
        group.max <- length(group.nums)
        group.nums <- c(group.nums, group.max + 1, group.max + 2)
        new.order <- order(match(memb, group.nums))
        legend.title <- 'Communities (#)'
      } else if (v.attr == 'lobe') {
        group.nums <- c(eval(parse(text=graph$atlas))[, levels(lobe)])
        group.max <- length(group.nums)
        group.nums <- c(group.nums, 'Inter', '')
        new.order <- order(memb)
        legend.title <- 'Lobe'
      }
      dat <- dat[new.order, new.order]
      dat.m <- melt(dat)
      setDT(dat.m)
      dat.m[, memb1 := group.nums[vertex_attr(graph, v.attr, as.character(Var1))]]
      dat.m[, memb2 := group.nums[vertex_attr(graph, v.attr, as.character(Var2))]]
      dat.m[, memb := ifelse(value == 1,
                             ifelse(memb1 == memb2, memb1, group.nums[group.max + 1]),
                             group.nums[group.max + 2])]
      if (v.attr == 'comm') {
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
      ggtitle(group) +
      theme(axis.ticks=element_blank(),
           axis.text.x=element_text(size=0.7*base_size, angle=45,
                                    color=corrs.m[1:Nv, color.text]),
           axis.title.x=element_blank(),
           axis.text.y=element_text(size=0.7*base_size,
                                    color=corrs.m[rev(1:Nv), color.text]),
           axis.title.y=element_blank()) +
      labs(fill=corrs.m[, unique(legend.t)]) +
      ylim(rev(levels(corrs.m$Var2)))

  } else {
    corrs.m <- melt(corrs)
    matplot <- ggplot(corrs.m, aes(Var1, Var2, fill=value)) +
      geom_tile() +
      ggtitle(group) +
      theme(axis.ticks=element_blank(),
           axis.text.x=element_text(size=0.7*base_size, angle=45),
           axis.title.x=element_blank(),
           axis.text.y=element_text(size=0.7*base_size),
           axis.title.y=element_blank()) +
      ylim(rev(levels(corrs.m$Var2)))
    if (identical(sum(abs(corrs)) - sum(corrs == 1), 0)) {
      matplot <- matplot + scale_fill_gradient2(low='white', high='blue') +
        theme(legend.position='none')
    }
  }
  return(matplot)
}
