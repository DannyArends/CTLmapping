\name{QTLimage}
\alias{QTLimage}

\title{Plot QTL curves}

\description{
  Plots the QTL heatmap of a genome wide QTL scan (part of the output of \code{\link{CTLscan}}).
}

\usage{
QTLimage(x, do.grid=TRUE, grid.col = "black", verbose = FALSE, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{do.grid}{ Adds a grid to the plot. }
  \item{grid.col}{ Color of the grid lines. }
  \item{verbose}{ Be verbose. }
  \item{\dots}{ Passed to the function \code{\link[graphics]{plot}} when it is called. }
}
\value{
  None.
}

\details{
  None.
}

\examples{
  library(ctl)          # Load CTL library
  data(ath.result)      # Arabidopsis Thaliana results
  QTLimage(ath.result)  # Plot only the QTLimage
}

\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Scan for CTL
    \item \code{\link{print.CTLscan}} - Print a summary of a CTLscan
    \item \code{\link[graphics]{par}} - Plot parameters
    \item \code{\link[grDevices]{colors}} - Colors used in plotting
  }
}

\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\keyword{hplot}
