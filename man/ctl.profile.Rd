\name{CTLprofile}
\alias{CTLprofile}
\title{CTLprofile - Extract single CTL interaction profiles }

\description{
  Extract for this CTL interaction profile the corresponding: phenotype x marker (p2m matrix) 
}

\usage{
  CTLprofile(CTLobject, pheno.col=1, significance = 0.05, verbose = TRUE)
}

\arguments{
  \item{CTLobject}{ A single element from an object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{pheno.col}{ Phenotype from the CTLobject to use. }
  \item{verbose}{ Be verbose. }
}

\details{
  CTL thresholds can be obtained using the \code{\link{CTLpermute}}  
  
  These matrices can be combined with QTL information to perform de novo 
  reconstruction of interaction networks.
}

\value{
  Matrix: phenotypes x marker.
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
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  ctl_result <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes,1:3,n.perm=10)
  p2m_matrix <- CTLprofile(ctl_result)
}

\keyword{methods}
