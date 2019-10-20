#' Perform an analysis with random graphs for brain MRI data
#'
#' \code{analysis_random_graphs} performs the steps needed for doing typical
#' graph theory analyses with brain MRI data if you need to generate equivalent
#' random graphs. This includes calculating \emph{small world} parameters and
#' normalized \emph{rich club} coefficients.
#'
#' \code{analysis_random_graphs} does the following:
#' \enumerate{
#'   \item Generate \code{N} random graphs for each graph and density/threshold
#'   \item Write graphs to disk in \code{savedir}. Read them back into \code{R}
#'     and combine into lists; then write these lists to disk. You can later
#'     delete the individual \code{.rds} files afterwards.
#'   \item Calculate \emph{small world} parameters, along with values for a few
#'     global graph measures that may be of interest.
#'   \item Calculate \emph{normalized rich club coefficients} and associated
#'     p-values.
#' }
#'
#' @param g.list List of \code{brainGraphList} objects; the length of this list
#'   should equal the number of thresholds/densities in the study
#' @param savedir Character string specifying the directory in which to save the
#'   generated graphs. Default: current working directory
#' @inheritParams Creating_Graphs
#' @export
#'
#' @return \code{analysis_random_graphs} returns a \emph{list} containing:
#' \item{rich}{A data table containing normalized rich-club coefficients and
#'   p-values}
#' \item{small}{A data table with small-world parameters}
#' \item{rand}{A data table with some global graph measures for all random
#'   graphs generated}
#'
#' @name Random Graphs
#' @rdname random_graphs
#' @family Random graph functions
#' @seealso \code{\link{small.world}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' rand_all <- analysis_random_graphs(g.norm, 1e2,
#'   savedir='/home/cwatson/dti/rand', clustering=F)
#' }

analysis_random_graphs <- function(g.list, level=g.list[[1]]$level, N=100L, savedir='.', ...) {
  Group <- NULL

  # Check if components are 'brainGraphList' objects
  matches <- vapply(g.list, inherits, logical(1), 'brainGraphList')
  if (any(!matches)) stop("Input must be a list of 'brainGraphList' objects.")

  kNumThresh <- length(g.list)
  kNumSubj <- vapply(g.list, function(x) length(x$graphs), integer(1))
  if (!dir.exists(savedir)) dir.create(savedir, recursive=TRUE)

  # Loop through all graphs, generate random graphs, calculate rich club coeff's
  phi.norm <- rand.dt <- small.dt <- vector('list', kNumThresh)
  for (i in seq_along(g.list)) {
    phi.norm[[i]] <- vector('list', kNumSubj[i])

    print(paste0('Random graphs for threshold #', i, '; ', format(Sys.time(), '%H:%M:%S')))
    progbar <- txtProgressBar(min=0, max=kNumSubj[i], style=3)
    for (j in seq_along(g.list[[i]]$graphs)) {
      rand <- sim.rand.graph.par(g.list[[i]][j], level=level, N=N, atlas=g.list[[i]]$atlas, ...)
      saveRDS(rand, file=file.path(savedir, sprintf('rand%02i_sub%04i%s', i, j, '.rds')))
      phi.norm[[i]][[j]] <- rich_club_norm(g.list[[i]][j], rand=rand)
      rm(list='rand')
      gc()
      setTxtProgressBar(progbar, j)
    }
  }
  close(progbar)

  # Get all small-worldness and random measures, all thresholds
  #-----------------------------------------------------------------------------
  for (i in seq_along(g.list)) {
    fnames <- list.files(savedir, sprintf('rand%02i.*', i), full.names=TRUE)
    rand.all <- lapply(fnames, readRDS)
    rand.all <- as_brainGraphList(rand.all, type='random', level=level)
    saveRDS(rand.all, file=file.path(savedir, sprintf('rand%02i%s', i, '_all.rds')))

    small.dt[[i]] <- small.world(g.list[[i]], rand.all)
    rand.dt[[i]] <- get_rand_attrs(rand.all, level)
    rm(list='rand.all')
    gc()
  }
  rand.dt <- rbindlist(rand.dt)
  small.dt <- rbindlist(small.dt)
  setkey(small.dt, Group, density)
  rich.dt <- rbindlist(lapply(phi.norm, rbindlist))

  return(list(rich=rich.dt, small=small.dt, rand=rand.dt))
}

#' Convenience function to get attributes for lists of random graphs
#'
#' @param bg.list A \code{brainGraphList} object, created in
#'   \code{\link{analysis_random_graphs}}
#' @keywords internal
#' @return A \code{data.table}

get_rand_attrs <- function(bg.list, level) {
  threshold <- NULL
  rand <- bg.list$graphs

  grps <- groups(bg.list)
  N <- lengths(rand)

  mod <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'mod'))
  Cp <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'Cp'))
  Lp <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'Lp'))
  E.global <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'E.global'))

  if (level == 'subject') {
    ids <- names(rand)
    DT <- data.table(Study.ID=rep(ids, times=N), Group=rep(grps, times=N))
  } else if (level == 'group') {
    DT <- data.table(Group=rep(grps, times=N))
  }
  DT <- cbind(DT, data.table(mod=mod, Cp=Cp, Lp=Lp, E.global=E.global))
  if ('density' %in% graph_attr_names(rand[[1]][[1]])) {
    densities <- vapply(rand, function(x) graph_attr(x[[1]], 'density'), numeric(1))
    DT[, density := rep(densities, times=N)]
  }
  if ('threshold' %in% graph_attr_names(rand[[1]][[1]])) {
    threshes <- vapply(rand, function(x) graph_attr(x[[1]], 'threshold'), numeric(1))
    DT[, threshold := rep(threshes, times=N)]
  }
  return(DT)
}
