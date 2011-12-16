\name{plot.QCLpermute}
\alias{plot.QCLpermute}

\title{Differential correlation versus Likelihood curves}

\description{
  Differential correlation versus Likelihood curves (the output of \code{\link{QCLpermute}}).
}

\usage{
\method{plot}{QCLpermute}(x, \dots)
}			    

\arguments{
  \item{x}{ An object of class \code{"QCLpermute"}, as output by \code{\link{QCLpermute}}. }
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
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes,1, n.perm=100)
  plot(qcl_result[[1]]$p)
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
