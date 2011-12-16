\name{plot.QCLobject}
\alias{plot.QCLobject}

\title{Plot QCL curves or heatmaps}

\description{
  Plot the QCL curve or heatmaps for a genome scan (the output of \code{\link{QCLscan}}).
}

\usage{
\method{plot}{QCLobject}(x, \dots)
}			    

\arguments{
  \item{x}{ An object of class \code{"QCLobject"}, as output by \code{\link{QCLscan}}. }
  \item{\dots}{ Passed to the function \code{\link{image.QCLobject}} \code{\link{plot.QCLscan}} when it is called. }
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
  plot(qcl_result)
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
