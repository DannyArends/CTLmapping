\name{CTLnetwork}
\alias{CTLnetwork}
\title{CTLnetwork - Interaction network from a genome-wide CTLscan of multiple traits }
\description{
  Create a file containing the interaction network from a genome-wide CTLscan of multiple traits.
}
\usage{
CTLnetwork(CTLobject, mapinfo, significance=0.05, LODdrop=2, what=c("names","ids"), add.qtls=FALSE, file="")
}
\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{ The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the 
    location of the marker in centiMorgans and the 3rd column "Mbp" - The location of the marker 
    in Mega basepairs. If supplied the marker names (rownames) should match those in the 
    CTLobject (only significant markers will be annotated). }
  \item{significance}{ Significance threshold for a genome wide false discovery rate (FDR). }
  \item{LODdrop}{ Drop in LOD score needed before we assign an edge type. }
  \item{what}{ Return trait and marker names or column numbers (for indexing). }
  \item{add.qtls}{ Should marker QTL trait interactions be added to the generated sif network 
    file, QTLs are included when they are above -log10(significance/n.markers). }
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
  A matrix with significant CTL interactions and information in 5 Columns:
  \itemize{
    \item TRAIT1 - Trait ID of the origin trait
    \item MARKER - Marker ID at which the CTL was found
    \item TRAIT2 - Trait ID of the target trait
    \item LOD_C  - LOD score of the CTL interaction
    \item CAUSAL - Type of edge determined by QTL LOD-drop:
      \itemize{
        \item NA - CTL/QTL for TRAIT1 and/or TRAIT2 not found
        \item -1 - TRAIT1 is DOWNSTREAM of TRAIT2
        \item 0  - UNDETERMINED Edge
        \item 1  - TRAIT1 is UPSTREAM of TRAIT2
      }
    \item LOD_T1 - QTL LOD-score of TRAIT1 at MARKER
    \item LOD_T2 - QTL LOD-score of TRAIT2 at MARKER
  }
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
  data(ath.result)       # Arabidopsis Thaliana results
  data(ath.metabolites)  # Arabidopsis Thaliana data set

  ctls <- CTLnetwork(ath.result, significance = 0.1)
  op <- par(mfrow = c(2,1))
  plot(ctls)
  ctl.lineplot(ath.result, ath.metab$map, significance=0.1)
}
\keyword{methods}
