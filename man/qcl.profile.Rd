\name{QCLprofile}
\alias{QCLprofile}
\title{QCLprofile - Extract single QCL interaction profiles }

\description{
	Extract for this QCL interaction profile the corresponding: phenotype x marker (p2m matrix) 
}

\usage{
  QCLprofile(QCL, qcl.threshold=0.4, verbose=TRUE)
}

\arguments{
  \item{QCL}{ A single element from an object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant. (See \code{\link{QCLpermute}}) }
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
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes)
  p2m_matrix <- QCLprofile(qcl_result[[2]])
}

\keyword{methods}
