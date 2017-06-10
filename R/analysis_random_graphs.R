#' Perform an analysis with random graphs for brain MRI data
#'
#' This function is not quite a "proper" function. It performs the steps needed
#' for doing typical graph theory analyses with brain MRI data if you need to
#' generate equivalent random graphs. This includes calculating \emph{small
#' world} parameters and normalized \emph{rich club} coefficients.
#'
#' The steps that are performed are:
#' \enumerate{
#'   \item \code{N} random graphs are generated for each group and
#'   density/threshold (and \emph{subject} if you have subject-specific graphs).
#'   \item These graphs are all written to disk in \code{savedir}. All of these
#'   are read back into \code{R} and combined into lists; these lists are also
#'   written to disk (in a sub-directory named \code{ALL}), so you can delete
#'   the individual \code{.rds} files afterwards.
#'   \item \emph{Small world} parameters are calculated, along with values for
#'   a few global graph measures that may be of interest.
#'   \item \emph{Normalized rich club coefficients} and associated p-values will
#'   be calculated.
#' }
#'
#' @param g.list List of lists containing \code{igraph} graph objects
#' @param N Integer specifying number of random graphs to generate per
#'   individual graph (default: 100)
#' @param savedir Character string specifying the directory in which to save the
#'   generated graphs (default: current working directory)
#' @param ... Other arguments passed to \code{\link{sim.rand.graph.par}} (e.g.
#'   \code{clustering=T})
#' @export
#'
#' @return A list containing:
#' \item{rich}{A list object containing normalized rich-club coefficients and
#'   p-values}
#' \item{small}{A data table with small-world parameters}
#' \item{rand}{A data table with some global graph measures for all random
#'   graphs generated}
#'
#' @family Null graph functions
#' @seealso \code{\link{small.world}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' rand_all <- analysis_random_graphs(g.norm, 1e2,
#'   savedir='/home/cwatson/dti/rand', clustering=F)
#' }

analysis_random_graphs <- function(g.list, N=100, savedir='.', ...) {
  Group <- Study.ID <- threshold <- NULL
  if (!isTRUE(dir.exists(paste0(savedir, '/ALL')))) {
    dir.create(paste0(savedir, '/ALL/'), recursive=TRUE)
  }

  # Convenience function
  create.rand.dt <- function(rand, group, V1, N) {
    rand.mod <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'mod'))
    rand.Cp <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'Cp'))
    rand.E.glob <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'E.global'))
    rand.Lp <- unlist(lapply(rand, vapply, graph_attr, numeric(1), USE.NAMES=FALSE, 'Lp'))
    rand.dt <- data.table(V1=rep(V1, each=N),
                          Group=rep(group, N), mod=rand.mod,
                          Cp=rand.Cp, Lp=rand.Lp,
                          E.global=rand.E.glob)
    return(rand.dt)
  }
  nesting <- 2
  if (is_igraph(g.list[[1]][[1]][[1]])) nesting <- 3

  # Loop through all graphs, generate random graphs, calculate rich club coeff's
  phi.norm <- vector('list', length=length(g.list))
  groups <- vector('character', length=length(g.list))
  for (i in seq_along(g.list)) {
    phi.norm[[i]] <- vector('list', length=length(g.list[[i]]))

    if (nesting == 2) {
      groups[i] <- g.list[[i]][[1]]$Group
      print(paste0('Random graphs for group #', i, '; ',
                   format(Sys.time(), '%H:%M:%S')))
      progbar <- txtProgressBar(min=0, max=length(g.list[[i]]), style=3)
      for (j in seq_along(g.list[[i]])) {
        rand <- sim.rand.graph.par(g.list[[i]][[j]], N, ...)
        saveRDS(rand, file=paste0(savedir, '/',
                                  sprintf('rand%i_thr%02i%s', i, j, '.rds')))
        phi.norm[[i]][[j]] <- rich_club_norm(g.list[[i]][[j]], rand=rand)
        rm(list='rand')
        gc()
        setTxtProgressBar(progbar, j)
      }
    } else {
      groups[i] <- g.list[[i]][[1]][[1]]$Group
      kNumSubjs <- length(g.list[[i]][[1]])
      for (j in seq_along(g.list[[i]])) {
        phi.norm[[i]][[j]] <- vector('list', length=kNumSubjs)
        print(paste0('Random graphs for group #', i, '; threshold #', j, '; ',
                     format(Sys.time(), '%H:%M:%S')))
        progbar <- txtProgressBar(min=0, max=length(g.list[[i]][[j]]), style=3)

        for (k in seq_along(g.list[[i]][[j]])) {
          rand <- sim.rand.graph.par(g.list[[i]][[j]][[k]], N, ...)
          saveRDS(rand, file=paste0(savedir, '/',
                                    sprintf('rand%i_thr%02i_subj%03i%s', i, j, k, '.rds')))
          phi.norm[[i]][[j]][[k]] <- rich_club_norm(g.list[[i]][[j]][[k]], rand=rand[[k]])
          rm(list='rand')
          gc()
          setTxtProgressBar(progbar, k)
        }
      }
    }
  }

  # Get all small-worldness and random measures, all thresholds
  #-----------------------------------------------------------------------------
  rand.dt <- small.dt <- ids <- vector('list', length=length(g.list))
  if (nesting == 2) {
    densities <- sapply(g.list[[1]], graph_attr, 'density')
    for (i in seq_along(g.list)) {
      fnames <- list.files(savedir, sprintf('rand%i_thr.*', i), full.names=TRUE)
      rand.all <- lapply(fnames, readRDS)
      saveRDS(rand.all, file=paste0(savedir, '/ALL/rand', i, '_all.rds'))
      small.dt[[i]] <- small.world(g.list[[i]], rand.all)

      rand.dt[[i]] <- create.rand.dt(rand.all, groups[i], V1=densities, N)
      setnames(rand.dt[[i]], 1, 'density')
      setkey(rand.dt[[i]], density)
      rm(list='rand.all')
      gc()
    }
    rand.dt <- rbindlist(rand.dt)
    small.dt <- rbindlist(small.dt)
    small.dt[, Group := rep(groups, times=lengths(g.list))]
    setkey(small.dt, Group, density)

  } else {
    for (i in seq_along(g.list)) {
      if ('name' %in% graph_attr_names(g.list[[i]][[1]][[1]])) {
        ids[[i]] <- sapply(g.list[[i]][[1]], graph_attr, 'name')
      } else {
        ids[[i]] <- seq_along(g.list[[i]][[1]])
      }
      rand.dt[[i]] <- small.dt[[i]] <- vector('list', length=length(g.list[[i]]))
      for (j in seq_along(g.list[[i]])) {
        fnames <- list.files(savedir, sprintf('rand%i_thr%02i.*', i, j), full.names=TRUE)
        rand.all <- lapply(fnames, readRDS)
        saveRDS(rand.all, file=paste0(savedir, '/ALL/',
                                      sprintf('rand%i_thr%02i%s', i, j, '.rds')))
        small.dt[[i]][[j]] <- small.world(g.list[[i]][[j]], rand.all)
        small.dt[[i]][[j]]$threshold <- j
        small.dt[[i]][[j]]$Study.ID <- ids[[i]]
        rand.dt[[i]][[j]] <- create.rand.dt(rand.all, groups[i], V1=ids[[i]], N)
        setnames(rand.dt[[i]][[j]], 1, 'Study.ID')
        rand.dt[[i]][[j]]$threshold <- j
        rm(list='rand.all')
        gc()
      }
    }
    rand.dt <- rbindlist(lapply(rand.dt, rbindlist))
    rand.dt[, Group := as.factor(Group)]
    rand.dt[, Study.ID := as.factor(Study.ID)]
    setcolorder(rand.dt, c('threshold', 'Group', 'Study.ID', names(rand.dt)[3:6]))
    setkey(rand.dt, threshold, Group, Study.ID)

    small.dt <- rbindlist(lapply(small.dt, rbindlist))
    kNumSubjs <- sapply(g.list, function(x) length(x[[1]]))
    small.dt[, Group := rep(groups, times=kNumSubjs*length(g.list[[1]]))]
    small.dt[, Group := as.factor(Group)]
    small.dt[, Study.ID := as.factor(Study.ID)]
    setkey(small.dt, Group, threshold)
  }

  return(list(rich=phi.norm, small=small.dt, rand=rand.dt))
}
