#' Perform an analysis with random graphs for brain MRI data
#'
#' This function is not quite a "proper" function. It performs the steps needed
#' for doing typical graph theory analyses with brain MRI data if you need to
#' generate equivalent random graphs. This includes calculating \emph{small
#' world} parameters and normalized \emph{rich club} coefficients.
#'
#' First, a number of random graphs are generated for each group and
#' density/threshold (and subject if you have subject-specific graphs). These
#' graphs are all written to disk in a location you specify. All of these are
#' read back into \code{R} and combined into large lists; these large lists are
#' also written to disk (in a sub-directory named \code{ALL}), so you can delete
#' the individual \code{.rds} files afterwards. Once all of the random graphs
#' have been generated, the \emph{small world} parameters are calculated, along
#' with values for a few global graph measures that may be of interest.
#' Additionally, the \emph{normalized rich club coefficients} and associated
#' p-values will be calculated.
#'
#' @param g A list object containing all graphs (may be nested)
#' @param N Integer specifying number of random graphs to generate per
#'   individual graph
#' @param covars Data table of covariates (used for Group and subject names)
#' @param savedir Character string specifying the directory in which to save the
#'   generated graphs (default: current working directory)
#' @param ... Other arguments passed to \code{\link{sim.rand.graph.par}} (e.g.
#'   \emph{clustering=F})
#' @export
#'
#' @return A list containing:
#' \item{rich}{A list object containing normalized rich-club coefficients and
#'   p-values}
#' \item{small}{A data table with small-world parameters}
#' \item{rand}{A data table with some global graph measures for all random
#'   graphs generated}
#'
#' @seealso \code{\link{sim.rand.graph.par}, \link{small.world},
#'   \link{rich.club.norm}}
#' @author Christopher G. Watson, \email{cgwatson@@bu.edu}
#' @examples
#' \dontrun{
#' rand_all <- random_graph_analysis(g.norm, 1e2, covars.dti,
#'   savedir='~/dti/rand', clustering=F)
#' }

analysis_random_graphs <- function(g, N, covars, savedir='.', ...) {
  Group <- Study.ID <- threshold <- NULL
  if (!isTRUE(dir.exists(paste0(savedir, '/ALL')))) {
    dir.create(paste0(savedir, '/ALL/'), recursive=TRUE)
  }

  # Convenience function
  create.rand.dt <- function(rand, group, V1, kNumRand) {
    rand.mod <- sapply(rand, sapply, graph_attr, 'mod')
    rand.E.global <- sapply(rand, sapply, graph_attr, 'E.global')
    rand.Cp <- sapply(rand, sapply, graph_attr, 'Cp')
    rand.Lp <- sapply(rand, sapply, graph_attr, 'Lp')
    rand.dt <- data.table(V1=rep(V1, each=kNumRand),
                          Group=rep(group, kNumRand), mod=as.vector(rand.mod),
                          Cp=as.vector(rand.Cp), Lp=as.vector(rand.Lp),
                          E.global=as.vector(rand.E.global))
    return(rand.dt)
  }
  # Check if the graphs have 3 or 2 nested levels
  nesting <- 2
  if (is.igraph(g[[1]][[1]][[1]])) nesting <- 3

  # Loop through all graphs, generate random graphs, calculate rich club coeff's
  phi.norm <- vector('list', length=length(g))
  for (i in seq_along(g)) {
    phi.norm[[i]] <- vector('list', length=length(g[[i]]))

    if (nesting == 2) {
      print(paste0('Random graphs for group #', i, '; ',
                   format(Sys.time(), '%H:%M:%S')))
      progbar <- txtProgressBar(min=0, max=length(g[[i]]), style=3)
      for (j in seq_along(g[[i]])) {
        rand <- sim.rand.graph.par(g[[i]][[j]], N, ...)
        saveRDS(rand, file=paste0(savedir, '/',
                                  sprintf('rand%i_thr%02i%s', i, j, '.rds')))
        phi.norm[[i]][[j]] <- rich.club.norm(g[[i]][[j]], rand=rand)
        rm(rand)
        gc()
        setTxtProgressBar(progbar, j)
      }
    } else {
      for (j in seq_along(g[[i]])) {
        phi.norm[[i]][[j]] <- vector('list', length=length(g[[i]][[j]]))
        print(paste0('Random graphs for group #', i, '; threshold #', j, '; ',
                     format(Sys.time(), '%H:%M:%S')))
        progbar <- txtProgressBar(min=0, max=length(g[[i]][[j]]), style=3)

        for (k in seq_along(g[[i]][[j]])) {
          rand <- sim.rand.graph.par(g[[i]][[j]][[k]], N, ...)
          saveRDS(rand, file=paste0(savedir, '/',
                                    sprintf('rand%i_thr%02i_subj%03i%s', i, j, k, '.rds')))
          phi.norm[[i]][[j]][[k]] <- rich.club.norm(g[[i]][[j]][[k]], rand=rand)
          rm(rand)
          gc()
          setTxtProgressBar(progbar, k)
        }
      }
    }
  }

  # Get all small-worldness and random measures, all thresholds
  #-----------------------------------------------------------------------------
  groups <- covars[, levels(Group)]
  fnames <- rand.dt <- small.dt <- vector('list', length=length(g))
  if (nesting == 2) {
    densities <- sapply(g[[1]], graph_attr, 'density')
    for (i in seq_along(g)) {
      fnames <- list.files(savedir, sprintf('rand%i_thr.*', i), full.names=T)
      rand_all <- lapply(fnames, readRDS)
      saveRDS(rand_all, file=paste0(savedir, '/ALL/rand', i, '_all.rds'))
      small.dt[[i]] <- small.world(g[[i]], rand_all)

      rand.dt[[i]] <- create.rand.dt(rand_all, groups[i], V1=densities, N)
      setnames(rand.dt[[i]], 1, 'density')
      setkey(rand.dt[[i]], density)
      rm(rand_all)
      gc()
    }
    rand.dt <- rbindlist(rand.dt)
    small.dt <- rbindlist(small.dt)
    small.dt[, Group := rep(groups, times=lengths(g))]
    setkey(small.dt, Group, density)

  } else {
    for (i in seq_along(g)) {
      fnames[[i]] <- rand.dt[[i]] <- small.dt[[i]] <- vector('list', length=length(g[[i]]))
      for (j in seq_along(g[[i]])) {
        fnames[[i]][[j]] <- list.files(savedir,
                                       sprintf('rand%i_thr%02i.*', i, j), full.names=T)
        rand.all <- lapply(fnames[[i]][[j]], readRDS)
        saveRDS(rand.all, file=paste0(savedir, '/ALL/',
                                      sprintf('rand%i_thr%02i%s', i, j, '.rds')))
        small.dt[[i]][[j]] <- small.world(g[[i]][[j]], rand.all)
        small.dt[[i]][[j]]$threshold <- j
        small.dt[[i]][[j]]$Study.ID <- covars[Group == groups[i], Study.ID]
        rand.dt[[i]][[j]] <- create.rand.dt(rand.all, groups[i],
                                            V1=covars[Group == groups[i], Study.ID], N)
        setnames(rand.dt[[i]][[j]], 1, 'Study.ID')
        rand.dt[[i]][[j]]$threshold <- j
        rm(rand.all)
        gc()
      }
    }
    rand.dt <- rbindlist(lapply(rand.dt, rbindlist))
    rand.dt[, threshold := as.factor(threshold)]
    rand.dt[, Group := as.factor(Group)]
    rand.dt[, Study.ID := as.factor(Study.ID)]
    setcolorder(rand.dt, c('threshold', 'Group', 'Study.ID', names(rand.dt)[3:6]))
    setkey(rand.dt, threshold, Group, Study.ID)

    small.dt <- rbindlist(lapply(small.dt, rbindlist))
    kNumSubjs <- sapply(g, function(x) length(x[[1]]))
    small.dt[, Group := rep(groups, times=kNumSubjs*length(g[[1]]))]
    small.dt[, Group := as.factor(Group)]
    small.dt[, threshold := as.factor(threshold)]
    small.dt[, Study.ID := as.factor(Study.ID)]
    setkey(small.dt, Group, threshold)
  }

  return(list(rich=phi.norm, small=small.dt, rand=rand.dt))
}
