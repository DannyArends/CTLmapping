\name{CTLprofiles}
\alias{CTLprofiles}
\title{CTLprofiles - Extract CTL interaction profiles }

\description{
  Extract the CTL interaction profiles: phenotype x marker (p2m matrix) 
  and phenotype x phenotype (p2p matrix) from a \code{\link{CTLscan}}.
}

\usage{
  CTLprofiles(CTLobject, against = c("markers","phenotypes"), significance = 0.05, verbose=FALSE, warn=TRUE)
}

\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{against}{ Plot the CTL against either: markers or phenotypes. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{verbose}{ Be verbose. }
  \item{warn}{ Warn when there are not enough permutations. }
}

\details{
  CTL thresholds can be obtained using the \code{\link{CTLpermute}}  
  
  These matrices can be combined with QTL information to perform de novo 
  reconstruction of interaction networks.
  
  The 'against' parameter is by default set to "markers" which returns a 
  phenotype x markers matrix (p2m matrix), which should be comparible to 
  the QTL profiles of the traits. 
  
  When the 'against' parameter is set to "phenotypes" a phenotype x phenotype
  matrix (p2p matrix) is returned, showing the interactions between the phenotypes.
}

\value{
  Matrix: phenotypes x marker or phenotypes x phenotypes
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
  p2m_matrix <- CTLprofiles(ctl_result, against="markers")
  p2p_matrix <- CTLprofiles(ctl_result, against="phenotypes")
}

\keyword{methods}
