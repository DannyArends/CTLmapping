\name{CTLnetwork}
\alias{CTLnetwork}
\title{CTLnetwork - Interaction network from a genome-wide CTLscan of multiple traits }
\description{
  Create a file containing the interaction network from a genome-wide CTLscan of multiple traits.
}
\usage{
CTLnetwork(CTLobject, significance=0.05, what=c("names","ids"), add.qtls=FALSE, mapinfo, file="")
}
\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{significance}{ Significance threshold for a genome wide false discovery rate (FDR). }
  \item{what}{ Return trait and marker names or column numbers (for indexing). }
  \item{add.qtls}{ Should marker QTL trait interactions be added to the generated sif network 
    file, QTLs are included when they are above -log10(significance/n.markers). }
  \item{mapinfo}{ The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the 
    location of the marker in centiMorgans and the 3rd column "Mbp" - The location of the marker 
    in Mega basepairs. If supplied the marker names (rownames) should match those in the 
    CTLobject (only significant markers will be annotated). }
  \item{file}{ A connection, or a character string naming the file to print to. If "" (the 
    default), CTLnetwork prints to the standard output connection, the console unless redirected 
    by sink. }
}
\details{
Outputs a sif network file, and a node attribute file:
\itemize{
  \item ctlnet<FILE>.sif - Shows CTL connections from Trait to Marker with edge descriptions
  \item ctlnet<FILE>.nodes - Attributes of the nodes (Traits and Genetic markers) nodes to this 
    file can be used to either color chromosomes, or add chromosome locations.
}
}
\value{
  NONE
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
  CTLnetwork(ath.result, 0.01)
}
\keyword{methods}
