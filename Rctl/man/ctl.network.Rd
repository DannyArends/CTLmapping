\name{CTLnetwork}
\alias{CTLnetwork}
\title{CTLnetwork - Interaction network from a genome-wide CTLscan of multiple traits }
\description{
  Create a file containing the interaction network from a genome-wide CTLscan of multiple traits.
}
\usage{
CTLnetwork(CTLobject, mapinfo, significance=0.05, add.qtl=TRUE, verbose=FALSE)
}

\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{mapinfo}{ The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the location of the marker in 
                  centiMorgans and the 3rd column "Mbp" - The location of the marker in Mega basepairs. }
  \item{add.qtl}{ Should QTL lines be added to the generated sif. }
  \item{verbose}{ Be verbose }
}
\details{
Output 2 sif network files, and a node attributes file:
\itemize{
  \item network_full.sif    - Shows all the significant CTL connections between traits
  \item network_summary.sif - Shows only a single summarized CTL connection between any 2 traits
  \item node_attributes.txt - Attributes of the genetic marker nodes to color chromosomes
}
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
