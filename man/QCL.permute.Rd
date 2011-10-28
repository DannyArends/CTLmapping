\name{QCLpermute}
\alias{QCLpermute}
\alias{permuteQCL}
\title{QCLpermute - Permutation routine for QCL mapping }
\description{
	Permutation routine for QCL mapping
}
\usage{
QCLpermute(genotypes, phenotypes, pheno.col, n.perm=10, directory="permutations", verbose=TRUE, \dots)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column should we use in permutation analysis }
  \item{n.perm}{ Number of permutations }
  \item{directory}{ Directory in which to save the output }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link{QCLscan}} function }
}
\details{
Permutation is one of the fundamental theorem a fundamental approach of 
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
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  qcl_perms <- QCLpermute(ath.metab$genotypes, ath.metab$phenotypes,1)
  qcl_perms
}

\keyword{methods}
