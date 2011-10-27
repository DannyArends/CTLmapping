\name{image.QCLscan}
\alias{image.QCLscan}

\title{Plot genome-wide QCL on multiple traits}

\description{
  Plot the QCLs for genome-wide QCL on multiple traits (the output of \code{\link{QCLscan}}).
}

\usage{
\method{image}{QCLscan}(x, qcl.threshold= 0.35, against = c("markers","phenotypes"), \dots)
}			    

\arguments{
  \item{x}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant. (See \code{\link{QCLpermute}}) }
  \item{against}{ Plot the QCL against either: markers or phenotypes. }
  \item{\dots}{ Passed to the function \code{\link[graphics]{image}} when it is called. }
}
\value{
  None.
}

\details{
  None.
}

\examples{
 #TODO
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
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\keyword{hplot}
