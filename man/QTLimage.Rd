\name{QTLimage}
\alias{QTLimage}

\title{Plot QTL curves}

\description{
  Plots the QTL heatmap of a genome wide QTL scan (part of the output of \code{\link{QCLscan}}).
}

\usage{
QTLimage(x, onlySignificant = FALSE, significance = 0.05, do.grid=TRUE, grid.col = "black", verbose = FALSE, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"QCLobject"}, as output by \code{\link{QCLscan}}. }
  \item{onlySignificant}{ Add only the significant QTL profiles to the image. }
  \item{significance}{ Minimum significance level for QTL profiles to be added. }
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
  library(qcl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes)
  QTLimage(qcl_result)
}

\seealso{
  \itemize{
    \item \code{\link{QCLscan}} - Scan for QCL
    \item \code{\link{print.QCLscan}} - Print a summary of a QCLscan
    \item \code{\link[graphics]{par}} - Plot parameters
    \item \code{\link[grDevices]{colors}} - Colors used in plotting
  }
}

\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\keyword{hplot}
