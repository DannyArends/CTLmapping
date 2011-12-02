\name{QCLpermute.cross}
\alias{QCLpermute.cross}
\title{QCLpermute.cross - Permutation routine for QCL mapping (R/qtl cross object) }
\description{
	Permutation routine for QCL mapping (using an R/qtl cross object)
}
\usage{
QCLpermute.cross(cross, pheno.col, n.perm=10, n.cores=2, directory="permutations", saveFiles=FALSE, verbose=FALSE, \dots)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype column should we use in permutation analysis }
  \item{n.perm}{ Number of permutations }
  \item{n.cores}{ Number of CPU cores to use during the analysis }
  \item{directory}{ Directory in which to save the output }
  \item{saveFiles}{ Save the permjutation matrices to HDD }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link{QCLscan}} function }
}
\details{
TODO
\itemize{
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system
}
}
\value{
QCLpermute object, a list holding QCLscan objects ontained at each permutation. To view significance thresholds use 
the \code{\link{print.QCLpermute}} function.
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
  \dontrun{
    require(qtl)
    data(multitrait) # Arabidopsis Thaliana (R/qtl cross object)
    qcl_perms <- QCLpermute.cross(multitrait,1)
    qcl_perms
  }
}

\keyword{methods}
