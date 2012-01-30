\name{qcl-package}
\alias{qcl}
\alias{QCL}
\docType{package}
\title{
  QCL - QCL mapping in experimental crosses
}
\description{
  Analysis of experimental crosses to identify genetic markers 
  associated with correlation changes in quantitative trait (QCL)
  The additional correlation information obtained can be combined 
  with QTL information to perform de novo reconstruction of 
  interaction networks.

  For more background information about the method we refer to the 
  methodology article published in XX (201X) :

  The R package is a basic iomplementation and it includes the 
  following core functionality:
  \itemize{
    \item \code{\link{QCLscan}} - Scan for QCL
    \item \code{\link{QCLpermute}} - Permutation routine for QCL mapping
    \item \code{\link{QCLnetwork}} - Create a QCL network from a \code{\link{QCLscan}}
    \item \code{\link{QCLprofiles}} - Extract QCL interaction profiles 
    \item \code{\link{image.QCLobject}} - Heatmap overview of a QCLscan
    \item \code{\link{plot.QCLscan}} - Plot the QCL curve for a single trait
  }
  
  Also provide are demos
}
\details{
  More detailed information and/or examples are given per function as needed.
  Some additional functionality:
  \itemize{
    \item \code{\link{QCLscan.cross}} - Use an R/qtl cross object with \code{\link{QCLscan}}
    \item \code{\link{QCLpermute.cross}} - Use an R/qtl cross object with \code{\link{QCLpermute}}
  }
}
\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com}\cr
  Contributions from: Bruno Tesson, Pjotr Prins and Ritsert C. Jansen
}
\references{
  \itemize{
    \item TODO
  }
}
\keyword{package}

\seealso{
  \itemize{
    \item \code{\link{QCLscan}} - Scan for QCL
    \item \code{\link{QCLpermute}} - Permutation routine for QCL mapping
  }
}
