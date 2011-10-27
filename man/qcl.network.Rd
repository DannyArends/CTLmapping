\name{QCLnetwork}
\alias{QCLnetwork}
\title{QCLnetwork - Network from a genome-wide QCLscan of multiple traits }
\description{
	Creates a sif network from a genome-wide QCLscan of multiple traits
}
\usage{
QCLnetwork(QCLscan, qcl.threshold=0.45, verbose = FALSE)
}

\arguments{
  \item{QCLscan}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant. (See \code{\link{QCLpermute}}) }
  \item{verbose}{ Be verbose }
}
\details{
TODO
}
\value{
 Vector of phenotype names inside the network
}
\references{ 
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  library(qcl)
  data(multitrait)
  multitrait <- fill.geno(multitrait)
  qcl_result <- QCLscan(multitrait)
  in_network <- QCLnetwork(qcl_result)
  in_network
}

\keyword{methods}
