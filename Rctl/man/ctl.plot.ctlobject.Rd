\name{plot.CTLobject}
\alias{plot.CTLobject}

\title{Plot CTL curves or heatmaps}

\description{
  Plot the CTL curve or heatmaps for a genome scan (the output of \code{\link{CTLscan}}).
}

\usage{
\method{plot}{CTLobject}(x, pheno.col=1:length(x), \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{pheno.col}{ Which phenotype column(s) should we plot. Defaults to creating an image of all phenotype columns }
  \item{\dots}{ Passed to the function \code{\link{image.CTLobject}} \code{\link{plot.CTLscan}} when it is called. }
}

\value{
  None.
}

\details{
  None.
}

\examples{
  library(ctl)
  data(ath.result) # Arabidopsis Thaliana dataset
  plot(ath.result)
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
