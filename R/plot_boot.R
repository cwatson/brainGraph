#' Plot global graph measures with shaded regions calculated from bootstrapping
#'
#' This function takes a list of \code{\link[boot]{boot}} objects (the number of
#' elements equals the number of subject groups) and plots the observed value
#' across all graph densities. It returns a list containing: a
#' \code{\link{data.table}} with standard errors and 95\% confidence intervals
#' at each density, and 2 \code{\link[ggplot2]{ggplot}} objects with shaded
#' regions surrounding the observed values.
#'
#' The 95\% confidence intervals are calculated using the normal approximation.
#'
#' @param boot.dt A \code{data.table} output from \code{\link{boot_global}}
#' @param ylabel A character string to place on the y-axis label (default: NULL)
#' @param alpha A numeric indicating the opacity for
#' \code{\link[ggplot2]{geom_ribbon}}
#' @param ... Other parameters passed to \code{\link[ggplot2]{geom_ribbon}}
#' @export
#'
#' @return A list with the following elements:
#' \item{p1}{A ggplot object with ribbon representing standard error}
#' \item{p2}{A ggplot object with ribbon representing 95\% confidence interval}
#' @seealso \code{\link{boot_global}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' boot.mod.plots <- plot_boot(boot.mod$dt, ylab='Modularity')
#' }

plot_boot <- function(boot.dt, ylabel=NULL, alpha=0.4, ...) {
  meas <- Group <- se <- ci.low <- ci.high <- NULL

  # First plot, using the SEM
  bootplot <- ggplot(boot.dt, aes(x=density, y=meas, col=Group)) +
    geom_line() +
    geom_ribbon(aes(ymin=meas-se, ymax=meas+se, fill=Group), alpha=alpha, ...) +
    ylab(ylabel)

  # Use the estimated normal 95% CI instead of se
  bootplot.ci <- ggplot(boot.dt, aes(x=density, y=meas, col=Group)) +
    geom_line() +
    geom_ribbon(aes(ymin=ci.low, ymax=ci.high, fill=Group), alpha=alpha, ...) +
    ylab(ylabel)

  return(list(p1=bootplot, p2=bootplot.ci))
}
