## ---- SimData ----
## Simulate a tractography data set, counts only
## I'm going for the general case, where all
## the imaging is done on a per subject basis,
## with merging into groups to happen later.
## Subjects are in their own folder
library(brainGraph)
library(tidyverse)
simdir <- "Sims"
dir.create(simdir, showWarnings = FALSE)

sd <- set.seed(23)
n <- 45
data("dk.scgm")

IDNums <- sample(1:200, size=n)

fnames <- file.path(simdir, paste0("Sub", formatC(IDNums, width=3, flag=0)), "counts.txt")

mkMat <- function(file, regions, min=0, max=1000) {
  mat <- matrix(runif(regions*regions, min=min, max=1000), nrow=regions)
  mat <- pmax(mat, t(mat))
  diag(mat) <- 0
  dir.create(dirname(file), showWarnings = FALSE)
  write.table(mat, file, row.names = FALSE, col.names = FALSE)
}

mats <- lapply(fnames, mkMat, regions=nrow(dk.scgm))

set.seed(sd)

## ---- StudyDesign ----

countfiles <- list.files(path=simdir, pattern="counts.txt", full.names=TRUE, recursive=TRUE)

studydata <- tibble(countfile=countfiles)
studydata <- mutate(studydata,
                    Index=1:nrow(studydata),
                    ID=basename(dirname(countfile)),
                    IDNum=as.numeric(stringr::str_remove(ID, "Sub")),
                    group=factor(rep(c("A", "B", "C"), n/3)))
studydata

## ---- LoadData3Threshold ----

thresholds <- c(0.1, 0.2, 0.3)
VersionA <- create_mats(studydata$countfile, modality = "dti", threshold.by = "density", mat.thresh=thresholds)

## ---- LoadData3ThresholdStructure ----
names(VersionA)
dim(VersionA$A)
dim(VersionA$A.norm)
length(VersionA$A.norm.sub)

dim(VersionA$A.norm.sub[[1]])

## ---- LoadDataConsensus3Group ----
groupindexes <- list()
groupindexes[[1]] <- which(studydata$group=="A")
groupindexes[[2]] <- which(studydata$group=="B")
groupindexes[[3]] <- which(studydata$group=="C")
sub.thresh <- 0.25
VersionB <- create_mats(studydata$countfile, modality = "dti", threshold.by = "consensus", sub.thresh=0.25, inds=groupindexes, mat.thresh=thresholds)

## ---- LoadDataConsensus3GroupStructure ----
names(VersionB)
dim(VersionB$A)
dim(VersionB$A.norm)
length(VersionB$A.norm.sub)

dim(VersionB$A.norm.sub[[1]])
range(VersionB$A.norm.sub[[1]] -- VersionA$A.norm.sub[[1]])

## ---- CreateBrainGraphParams ----

params <- as.tibble(expand.grid(idx=studydata$Index, threshnum=1:length(thresholds)))
params <- mutate(params, 
                 thresh=thresholds[threshnum],
                 ID=studydata$ID[idx],
                 Group=studydata$group[idx])
params
## ---- CreateBrainGraph ----

# Wrapper function that can extract relevant data from 
# the params table
doOne <- function(pidx, params, Mats) {
  ## get the parameters from the table - force drop for tibbles
  idx <- params[pidx, "idx", drop=TRUE]
  thr <- params[pidx, "thresh", drop=TRUE]
  threshnum <- params[pidx, "threshnum", drop=TRUE]
  ID <- params[pidx, "ID", drop=TRUE]
  gr <- graph_from_adjacency_matrix(Mats$A.norm.sub[[threshnum]][,,idx],
                                    mode='undirected', diag=F,
                                    weighted=T)
  g.tmp <- set_brainGraph_attr(gr, "dk.scgm", modality='dti',
                               weighting='counts', threshold=thr,
                               subject=ID,
                               use.parallel=FALSE, A=Mats$A.norm.sub[[threshnum]][, , idx])
  return(g.tmp)
}


bgraphs <- parallel::mclapply(1:nrow(params), doOne, params=params, Mats=VersionA, mc.cores=4)
graph_attr_dt(bgraphs)
## ---- RandomGraphs ----

doOneRand <- function(gr, N=100, ...) {
  rand <- sim.rand.graph.par(gr, 
                             N, ...)
  phi.norm <- rich_club_norm(gr,  rand = rand)
  dns <- graph_attr(gr, "density")
  sw <- small.world(gr, rand)
  return(list(sw, phi.norm))
}

smw <- parallel::mclapply(bgraphs, doOneRand, N=10, mc.cores=4)

smw[[1]]