\name{QCL.permute}
\alias{QCL.permute}
\title{QCL.permute - Permutation routine for QCL mapping }
\description{
	Permutation routine for QCL mapping
}
\usage{
QCL.permute(cross, n.perm=10, directory="permutations", verbose=TRUE, \dots)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{n.perm}{ Number of permutations }
  \item{directory}{ Directory in which to save the output }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link{QCLscan}} function }
}
\details{
TODO
\itemize{
  \item NOTE: QCLscan is optimized for multicore machines its uses snow to use two cores
}
}
\value{
QCLpermute object, a list holding QCLscan objects from each permutation. To view significance thresholds use 
the \code{\link{significance.QCL}} function.
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
	#TODO
}

\keyword{methods}
