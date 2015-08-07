#' Plot correlation matrices
#'
#' This function will plot correlation matrices in the form of a ``heatmap''. It
#' will work for to two groups, and you have the choice to plot the vertices in
#' an order based on either community or lobe membership.
#'
#' @param c1 The correlation matrix for group 1
#' @param c2 The correlation matrix for group 2
#' @param ordered A logical indicating whether or not to order vertices
#' (default:TRUE)
#' @param type Character string, either 'comm' or 'lobe'
#' @param g1 An igraph graph object; not required if \emph{ordered} is FALSE
#' @param g2 An igraph graph object; not required if \emph{ordered} is FALSE
#' @param groups A character vector of group names
#' @export
#'
#' @return A list of ggplot objects
#' @seealso \code{\link[ggplot2]{geom_tile}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' mat.plots <- plot_corr_mat(corrs.group1[[N]]$r.thresh,
#'                            corrs.group2[[N]]$r.thresh, g1=g1[[N]],
#'                            g2=g2[[N]], groups=groups)
#' }

plot_corr_mat <- function(c1, c2, ordered=TRUE, type=c('comm', 'lobe'), g1,
                          g2, groups=c('Group 1', 'Group 2')) {
  base_size <- ifelse(nrow(c1) > 90, 7.5, 9)

  if (isTRUE(ordered)) {
    if (is.null(g1) | is.null(g2)) {
      stop('You must provide graph objects for vertex ordering!')
    }
    Nv <- nrow(c1)
    cols <- group.cols

    type <- match.arg(type)
      create.dt <- function(dat, graph, v.attr) {
        memb <- vertex_attr(graph, v.attr)
        if (v.attr == 'comm') {
          tab <- table(memb)
          group.nums <- as.integer(names(tab))
          group.max <- length(group.nums)
          group.nums <- c(group.nums, group.max + 1, group.max + 2)
          new.order <- order(match(memb, group.nums))
          legend.title <- 'Communities (#)'
        } else if (v.attr == 'lobe') {
          atlas <- graph$atlas
          atlas.dt <- eval(parse(text=data(list=atlas)))
          group.nums <- c(atlas.dt[, levels(lobe)])
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
      c1.m <- create.dt(c1, g1, type)
      c2.m <- create.dt(c2, g2, type)

      mats <- Map(function(w, y) {
        w$Group <- y
        ggplot(w, aes(Var1, Var2, fill=memb)) +
        geom_tile() +
        scale_fill_manual(values=w[, levels(color)]) +
        ggtitle(y) +
        theme(axis.ticks=element_blank(),
             axis.text.x=element_text(size=0.7*base_size, angle=45,
                                      color=w[1:Nv, color.text]),
             axis.title.x=element_blank(),
             axis.text.y=element_text(size=0.7*base_size,
                                      color=w[rev(1:Nv), color.text]),
             axis.title.y=element_blank()) +
        labs(fill=w[, unique(legend.t)]) +
        ylim(rev(levels(w$Var2)))},
        list(c1.m, c2.m), as.list(groups))

    return(list(g1=mats[[1]], g2=mats[[2]]))
  }
}
