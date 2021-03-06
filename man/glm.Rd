% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brainGraph_GLM.R, R/glm_methods.R
\name{GLM}
\alias{GLM}
\alias{brainGraph_GLM}
\alias{print.bg_GLM}
\alias{summary.bg_GLM}
\alias{plot.bg_GLM}
\alias{[.bg_GLM}
\title{Fit General Linear Models at each vertex of a graph}
\usage{
brainGraph_GLM(g.list, covars, measure, contrasts, con.type = c("t",
  "f"), outcome = NULL, X = NULL, con.name = NULL,
  alternative = c("two.sided", "less", "greater"), alpha = 0.05,
  level = c("vertex", "graph"), permute = FALSE,
  perm.method = c("freedmanLane", "terBraak", "smith", "draperStoneman",
  "manly", "stillWhite"), part.method = c("beckmann", "guttman",
  "ridgway"), N = 5000, perms = NULL, long = FALSE, ...)

\method{print}{bg_GLM}(x, ...)

\method{summary}{bg_GLM}(object, p.sig = c("p", "p.fdr", "p.perm"),
  contrast = NULL, alpha = object$alpha, digits = max(3L,
  getOption("digits") - 2L), print.head = TRUE, ...)

\method{plot}{bg_GLM}(x, region = NULL, which = c(1L:3L, 5L),
  ids = TRUE, ...)

\method{[}{bg_GLM}(x, i, j)
}
\arguments{
\item{g.list}{A \code{brainGraphList} object}

\item{covars}{A \code{data.table} of covariates}

\item{measure}{Character string of the graph measure of interest}

\item{contrasts}{Numeric matrix (for T statistics) or list of matrices (for F
statistics) specifying the contrast(s) of interest; if only one contrast is
desired, you can supply a vector (for T statistics)}

\item{con.type}{Character string; either \code{'t'} or \code{'f'} (for t or
F-statistics). Default: \code{'t'}}

\item{outcome}{Character string specifying the name of the outcome variable,
if it differs from the graph metric (\code{measure})}

\item{X}{Numeric matrix, if you wish to supply your own design matrix.
Ignored if \code{outcome != measure}.}

\item{con.name}{Character vector of the contrast name(s); if \code{contrasts}
has row/list names, those will be used for reporting results}

\item{alternative}{Character string, whether to do a two- or one-sided test.
Default: \code{'two.sided'}}

\item{alpha}{Numeric; the significance level. Default: 0.05}

\item{level}{Character string; either \code{vertex} (default) or
\code{graph}}

\item{permute}{Logical indicating whether or not to permute group labels.
Default: \code{FALSE}}

\item{perm.method}{Character string indicating the permutation method.
Default: \code{'freedmanLane'}}

\item{part.method}{Character string; the method of partitioning the design
matrix into covariates of interest and nuisance. Default: \code{'beckmann'}}

\item{N}{Integer; number of permutations to create. Default: \code{5e3}}

\item{perms}{Matrix of permutations, if you would like to provide your own.
Default: \code{NULL}}

\item{long}{Logical indicating whether or not to return all permutation
results. Default: \code{FALSE}}

\item{...}{Arguments passed to \code{\link{brainGraph_GLM_design}}}

\item{object, x}{A \code{bg_GLM} object}

\item{p.sig}{Character string specifying which P-value to use for displaying
significant results (default: \code{p})}

\item{contrast}{Integer specifying the contrast to plot/summarize; defaults
to showing results for all contrasts}

\item{digits}{Integer specifying the number of digits to display for P-values}

\item{print.head}{Logical indicating whether or not to print only the first
and last 5 rows of the statistics tables (default: \code{TRUE})}

\item{region}{Character string specifying which region's results to
plot; only relevant if \code{level='vertex'}. Default: \code{NULL}}

\item{which}{Integer vector indicating which of the 6 plots to print to the
plot device. Default: \code{c(1:3, 5)}}

\item{ids}{Logical indicating whether to plot subject ID's for outliers.
Otherwise plots the integer index}

\item{i}{Integer/character vector; the observation number(s) or row names to
select or remove}

\item{j}{Integer/character vector; the design matrix column number(s) or
names to select or remove}
}
\value{
An object of class \code{bg_GLM} containing some input-specific
  variables (\code{level}, \code{outcome}, \code{measure}, \code{con.type},
  \code{contrasts}, \code{con.name}, \code{alt}, \code{alpha},
  \code{permute}, \code{perm.method}, \code{part.method}, \code{N}) in
  addition to:
  \item{DT.Xy}{A data table from which the design matrices are created and
    the outcome variable, for all regions.}
  \item{X}{A named numeric matrix or a 3D array of the design matrix.
    Rownames are subject IDs, column names are predictor variables, and
    dimnames along the 3rd dimension are region names (if applicable). This
    is a 3D array only if \code{outcome != measure} and \code{level ==
    'vertex'}.}
  \item{y}{A named numeric matrix of the outcome variable. Rownames are Study
    IDs and column names are regions. There will be multiple columns only if
    \code{outcome == measure} and \code{level == 'vertex'}.}
  \item{DT}{A data table with an entry for each vertex (region) containing
    statistics of interest}
  \item{removed.subs}{A named integer vector in which the names are subject
    ID's of those removed due to incomplete data (if any). The integers
    correspond to the row number in the input \code{covars} table.}
  \item{runX}{If \code{outcome != measure} and \code{level == 'vertex'}, this
    will be a character vector of the regions for which the design matrix is
    invertible. Otherwise, it is \code{NULL}.}
  \item{runY}{Character vector of the regions for which the outcome variable
    has 0 variability. For example, if \code{level='vertex'} and
    \code{measure='degree'}, some regions may be disconnected or have the
    same degree for all subjects.}
  \item{atlas}{Character string of the atlas used (guessed based on the
    vertex count).}
  \item{perm}{A list containing: \emph{null.dist} (the null distribution of
    maximum statistics), \emph{thresh} (the statistic value corresponding
    to the \eqn{100 \times (1 - \alpha)}th\% percentile of the null
    distribution)}

The \code{plot} method returns a \emph{list} of \code{ggplot} objects
  (if installed) or writes the plots to a PDF in the current directory named
  \code{bg_GLM_diagnostics.pdf}

A \code{bg_GLM} object with the specified row(s) selected or removed
  from both \code{X} and \code{y}, and column(s) selected/removed from
  \code{X}
}
\description{
\code{brainGraph_GLM} specifies and fits a General Linear Model (GLM) at each
vertex for a given vertex measure (e.g. \emph{degree}) or at the graph-level
(e.g., \emph{global efficiency}). Given a contrast matrix or list of
contrast(s), and contrast type (for t- or F-contrast(s), respectively) it
will calculate the associated statistic(s) for the given contrast(s).

The \code{summary} method prints the results, only for which
\eqn{p < \alpha}, where \code{alpha} comes from the \code{bg_GLM} object.
\dQuote{Simple} P-values are used by default, but you may change this to the
FDR-adjusted or permutation P-values via the function argument \code{p.sig}.
You may also choose to subset by \emph{contrast}.

The \code{plot} method plots the GLM diagnostics (similar to that of
\code{\link[stats]{plot.lm}}). There are a total of 6 possible plots,
specified by the \code{which} argument; the behavior is the same as in
\code{\link[stats]{plot.lm}}. Please see the help for that function.

The \code{[} method allows you to select observations (i.e., rows of \code{X}
and \code{y}) and independent variables (i.e., columns of \code{X}) from a
\code{bg_GLM} object.
}
\details{
The \code{measure} argument will be the graph- or vertex-level measure of
interest. Often, this will serve as the model's \emph{outcome} (or dependent,
or response) variable; i.e., the variable typically denoted by \emph{y} in
GLMs. In other cases, you may wish to choose some other variable as the
outcome; e.g., IQ, age, etc. Then you could test for a direct association
between the network measure and outcome of interest, or test for another
association while adjusting for the network metric. For these applications,
you must provide the variable name via the \code{outcome} argument. This is
analogous to \code{-evperdat} in FSL's PALM and to \code{--pvr} in
FreeSurfer.
}
\note{
The \code{[} method is used when calculating \emph{studentized
residuals} and other \dQuote{leave-one-out} diagnostics, and typically should
not be called directly by the user.
}
\section{Design matrix}{

The GLM's \emph{design matrix} will often be identical to the \emph{model
matrix} associated with \code{lm} objects (if \dQuote{dummy} coding, the
default, is used) and is created from the input \code{data.table} and
arguments passed to \code{\link{brainGraph_GLM_design}}. The first column
should have the name of \code{getOption('bg.subject_id')} and its values
must match the \emph{name} graph-level attribute of the input graphs. The
covariates table must be supplied even if you provide your own design matrix
\code{X}.  If \code{level='vertex'} and \code{outcome == measure}, there will
be a single design for all regions but a separate model for each region
(since the graph measure varies by region). If \code{level='vertex'} and
\code{outcome != measure}, there will be a separate design (and, therefore, a
separate model) for each region even though the outcome is the same in all
models.
}

\section{Contrasts and statistics}{

Either t- or F-contrasts can be calculated (specified by \code{con.type}).
Multiple t-contrasts can be specified by passing a multi-row \emph{matrix} to
\code{contrasts}. Multiple F-contrasts can be specified by passing a
\emph{list} of matrices; all matrices must have the same number of columns.
All F-contrasts are necessarily \emph{two-sided}; t-contrasts can be any
direction, but only one can be chosen per function call.
If you choose \code{con.type="f"}, the calculated effect size is represented
by the \code{ESS} (\dQuote{extra sum of squares}), the additional variance
explained for by the model parameters of interest (as determined by the
contrast matrix). The standard error for F-contrasts is the sum of squared
errors of the \emph{full model}.
}

\section{Non-parametric permutation tests}{

You can calculate permutations of the data to build a null distribution of
the maximum statistic which corrects for multiple testing. To account for
complex designs, the design matrix must be \emph{partitioned} into covariates
of interest and nuisance; the default method is the \emph{Beckmann} method.
The default permutation strategy is that of Freedman & Lane (1983), and is
the same as that in FSL's \emph{randomise}. See \code{\link{randomise}}.
}

\examples{
\dontrun{
conmat <- matrix(c(0, 0, 0, 1), nrow=1)
rownames(conmat) <- 'Control > Patient'

res.lm <- brainGraph_GLM(g[[6]], covars=covars.all[tract == 1],
  measure='strength', contrasts=conmat, alt='greater', permute=TRUE, long=TRUE)
}
\dontrun{
## Save objects and then to multipage PDF
lmPlots <- plot(x)
ggsave('lmPlots.pdf', lmPlots)

## Save all the GLM sub-objects from MTPC analysis
res.mtpc <- mtpc(...)
glmPlots <- lapply(res.mtpc$res.glm, plot, which=1:6)
ml <- marrangeGrob(glmPlots, nrow=1, ncol=1)
ggsave('glmPlots.pdf', ml, width=8.5, height=11)
}
}
\seealso{
\code{\link[stats]{plot.lm}}

Other GLM functions: \code{\link{GLM design}},
  \code{\link{GLM fits}}, \code{\link{mtpc}}

Other Group analysis functions: \code{\link{Bootstrapping}},
  \code{\link{Mediation}}, \code{\link{NBS}},
  \code{\link{brainGraph_permute}}, \code{\link{mtpc}}
}
\author{
Christopher G. Watson, \email{cgwatson@bu.edu}
}
\concept{GLM functions}
\concept{Group analysis functions}
