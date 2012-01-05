\name{QCLpermute}
\alias{QCLpermute}
\alias{permuteQCL}
\title{QCLpermute - Permutation routine for QCL mapping }
\description{
  Permutation routine for QCL mapping
}
\usage{
QCLpermute(genotypes, phenotypes, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose=FALSE, \dots)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column should we use in permutation analysis }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{n.perm}{ Number of permutations }
  \item{n.cores}{ Number of CPU cores to use during the analysis }
  \item{genotype.values}{ Two values specifying the encoding of the genotype matrix }
  \item{directory}{ Directory in which to save the output }
  \item{saveFiles}{ Save the permjutation matrices to HDD }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link{QCLscan}} function }
}
\details{
Permutation is one of the fundamental approach to determine probability when nothing is known about the underlying NULL-ditribution
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
