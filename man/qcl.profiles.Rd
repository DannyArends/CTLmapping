\name{QCLprofiles}
\alias{QCLprofiles}
\title{QCLprofiles - Extract QCL interaction profiles }

\description{
	Extract the QCL interaction profiles: phenotype x marker (p2m matrix) 
  and phenotype x phenotype (p2p matrix) from a \code{\link{QCLscan}}.
}

\usage{
  QCLprofiles(QCLscan, qcl.threshold=0.4, against = c("markers","phenotypes"))
}

\arguments{
  \item{QCLscan}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant. (See \code{\link{QCLpermute}}) }
  \item{against}{ Plot the QCL against either: markers or phenotypes. }
}

\details{
  QCL thresholds can be obtained using the \code{\link{QCLpermute}}  
  
  These matrices can be combined with QTL information to perform de novo 
  reconstruction of interaction networks.
  
  The 'against' parameter is by default set to "markers" which returns a 
  phenotype x markers matrix (p2m matrix), which should be comparible to 
  the QTL profiles of the traits. 
  
  When the 'against' parameter is set to "phenotypes" a phenotype x genetic 
  marker matrix (p2m matrix) is returned, showing the number of markers at 
  which interactions take occur.
}

\value{
  Matrix: phenotypes x phenotypes or phenotypes x marker.
}
\references{ 
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  library(qcl)
  data(multitrait)
  multitrait <- fill.geno(multitrait)
  qcl_result <- QCLscan(multitrait)
  p2m_matrix <- QCLprofiles(qcl_result, against="markers")
  p2p_matrix <- QCLprofiles(qcl_result, against="phenotypes")
}

\keyword{methods}
