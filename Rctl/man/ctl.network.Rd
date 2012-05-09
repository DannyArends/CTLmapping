\name{CTLnetwork}
\alias{CTLnetwork}
\title{CTLnetwork - Interaction network from a genome-wide CTLscan of multiple traits }
\description{
  Create a file containing the interaction network from a genome-wide CTLscan of multiple traits.
}
\usage{
CTLnetwork(CTLobject, lod.threshold=5, chr.edges, add.qtl = TRUE, verbose = FALSE)
}

\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{lod.threshold}{ LOD score above which a CTL is included. (See \code{\link{CTLpermute}}) }
  \item{chr.edges}{ LOD score above which a CTL is included. (See \code{\link{CTLpermute}}) }
  \item{add.qtl}{ Should QTL lines be added to the generated sif. }
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
