\name{ctl.load}
\alias{ctl.load}
\title{ctl.load - Load CTLs calculated by the D2.0 version }
\description{
  Load CTLs calculated by the D2.0 version
}
\usage{
  ctl.load(genotypes = "ngenotypes.txt", phenotypes = "nphenotypes.txt", output = "ctlout", from=1, to, 
           verbose = FALSE)
}
\arguments{
  \item{genotypes}{ Original datafile containing the genotypes scanned. }
  \item{phenotypes}{ Original datafile containing the phenotypes scanned. }
  \item{output}{ Directory containing the output files. }
  \item{from}{ Start loading at which phenotype. }
  \item{to}{ Continue loading untill this phenotype. }
  \item{verbose}{ Be verbose. }
}
\details{
  TODO
}
\value{
  CTLobject, a list with at each index a CTLscan object:
  \itemize{
    \item $ctls - Matrix of differential correlation scores for each trait at each marker
    \item $qtl  - Vector of QTL lodscores for each marker (if a QTL scan was perfomed -qtl)
    \item $p    - Vector of maximum scores per marker obtained during permutations
    \item $l    - Matrix of LOD scores for CTL likelihood
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
  library(ctl)         # Load CTL library
}
\keyword{methods}
