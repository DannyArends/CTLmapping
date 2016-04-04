\name{CTLregions}
\alias{CTLregions}
\title{CTLregions - Get all significant interactions from a genome-wide CTLscan }
\description{
  Get all significant interactions from a genome-wide CTLscan.
}
\usage{
CTLregions(CTLobject, map, phenocol = 1, significance = 0.05, verbose = TRUE)
}

\arguments{
  \item{CTLobject}{         An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{           The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the 
                            location of the marker in centiMorgans and the 3rd column "Mbp" - The location 
                            of the marker in Mega basepairs. If supplied the marker names (rownames) 
                            should match those in the CTLobject. }
  \item{phenocol}{          Which phenotype column should we analyse. }
  \item{significance}{      Significance threshold to set a genome wide False Discovery Rate (FDR). }
  \item{verbose}{           Be verbose. }
}
\details{
  TODO
}
\value{
  A matrix significant CTL interactions with 4 columns: trait, marker, trait, lod
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
  data(ath.metab)
  data(ath.result)
  regions <- CTLregions(ath.result, ath.metab$map)
}

\keyword{methods}
