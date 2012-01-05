\name{QCLnetwork}
\alias{QCLnetwork}
\title{QCLnetwork - Interaction network from a genome-wide QCLscan of multiple traits }
\description{
  Create a file containing the interaction network from a genome-wide QCLscan of multiple traits.
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
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  library(qcl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes, 1:3, n.perm=10)
  in_network <- QCLnetwork(qcl_result)
  in_network
}

\keyword{methods}
