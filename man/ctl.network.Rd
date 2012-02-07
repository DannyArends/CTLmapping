\name{CTLnetwork}
\alias{CTLnetwork}
\title{CTLnetwork - Interaction network from a genome-wide CTLscan of multiple traits }
\description{
  Create a file containing the interaction network from a genome-wide CTLscan of multiple traits.
}
\usage{
CTLnetwork(CTLscan, ctl.threshold=0.45, verbose = FALSE)
}

\arguments{
  \item{CTLscan}{ An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{ctl.threshold}{ Threshold at which an CTL is considered significant. (See \code{\link{CTLpermute}}) }
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
  library(ctl)
  data(ath.result)
  in_network <- CTLnetwork(ath.result)
  in_network
}

\keyword{methods}
