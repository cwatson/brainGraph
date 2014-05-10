\name{plot.over.brain}
\alias{plot.over.brain}
\title{Draw an axial outline of a brain.}
\usage{
\method{plot}{over.brain}(flag = 0, fname = NULL)
}
\arguments{
  \item{flag}{binary indicating whether or not a png file
  should be saved}

  \item{fname}{the name of the file to be saved}
}
\description{
This function draws an outline of a brain (axial view), to
plot the nodes of a graph over it. It will optionally write
to a filename for output. It requires that the
\code{\link{png}} package be loaded. This specific png is
412x480 in resolution; the outline goes from ~x=50 -> 350,
y=60 -> 420.
}

