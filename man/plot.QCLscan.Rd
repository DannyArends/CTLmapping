\name{plot.QCLscan}
\alias{plot.QCLscan}

\title{Plot QCL curves}

\description{
  Plot the QCL curve for a genome scan (the output of \code{\link{QCLscan}}).
}

\usage{
\method{plot}{QCLscan}(x, pheno.col=1, qcl.threshold=0.3, do.legend = TRUE, \dots)
}			    

\arguments{
  \item{x}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{pheno.col}{ Phenotype to plot. }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant. (See \code{\link{QCLpermute}}) }
  \item{do.legend}{ Adds a legend with the different qcl.significance thresholds used. }
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
