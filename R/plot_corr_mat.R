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

plot_corr_mat <- function(c1, c2, ordered=TRUE, type=c('comm', 'lobe'), g1=NULL,
                          g2=NULL, groups=c('Group 1', 'Group 2')) {
  base_size <- ifelse(nrow(c1) > 100, 7.5, 9)

  if (isTRUE(ordered)) {
    if (is.null(g1) | is.null(g2)) {
      stop('You must provide graph objects for vertex ordering!')
    }
    type <- match.arg(type)
    if (type == 'comm') {
      x1 <- V(g1)$comm
      comms1 <- table(x1)
      y1 <- as.integer(names(comms1[rev(order(comms1))]))
      ord1 <- order(match(x1, y1))
      x2 <- V(g2)$comm
      comms2 <- table(x2)
      y2 <- as.integer(names(comms2[rev(order(comms2))]))
      ord2 <- order(match(x2, y2))
    } else if (type == 'lobe') {
      ord1 <- ord2 <- order(V(g1)$lobe)
    }
    c1 <- c1[ord1, ord1]
    c2 <- c2[ord2, ord2]
  }
  mats <- Map(function(x, y) {
    dat <- melt(x)
    dat$Group <- y
    ggplot(melt(x), aes(Var1, Var2, fill=value)) +
    geom_tile() +
    scale_fill_gradient2(low='white', high='blue') +
    ggtitle(y) +
    theme(legend.position='none', axis.ticks=element_blank(),
         axis.text.x=element_text(size=0.7*base_size, angle=45),
         axis.title.x=element_blank(),
         axis.text.y=element_text(size=0.7*base_size),
         axis.title.y=element_blank()) +
    ylim(rev(levels(melt(x)$Var2)))},
    list(c1, c2), as.list(groups))
  return(list(g1=mats[[1]], g2=mats[[2]]))
}
