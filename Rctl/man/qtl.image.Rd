\name{qtlimage}
\alias{qtlimage}
\title{Plot a QTL heatmap of the phenotypes scanned by CTLscan}
\description{
  Plots the QTL heatmap of a genome wide QTL scan (part of the output of \code{\link{CTLscan}}).
}
\usage{
qtlimage(x, do.grid=TRUE, grid.col="white", verbose=FALSE, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{do.grid}{ When TRUE, grid lines are added to the plot. }
  \item{grid.col}{ Color used for the grid lines, only used when do.grid = TRUE. }
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
  qtlimage(ath.result)  # Plot only the qtls
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
