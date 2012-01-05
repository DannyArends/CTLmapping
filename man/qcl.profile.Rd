\name{QCLprofile}
\alias{QCLprofile}
\title{QCLprofile - Extract single QCL interaction profiles }

\description{
  Extract for this QCL interaction profile the corresponding: phenotype x marker (p2m matrix) 
}

\usage{
  QCLprofile(QCLobject, pheno.col=1, significance = 0.05, verbose = TRUE)
}

\arguments{
  \item{QCLobject}{ A single element from an object of class \code{"QCLobject"}, as output by \code{\link{QCLscan}}. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{pheno.col}{ Phenotype from the QCLobject to use. }
  \item{verbose}{ Be verbose. }
}

\details{
  QCL thresholds can be obtained using the \code{\link{QCLpermute}}  
  
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
  library(qcl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes,1:3,n.perm=10)
  p2m_matrix <- QCLprofile(qcl_result)
}

\keyword{methods}
