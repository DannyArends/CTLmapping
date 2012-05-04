\name{CTLpermute}
\alias{CTLpermute}
\alias{permuteCTL}
\title{CTLpermute - Permutation routine for CTL mapping }
\description{
  Permutation routine for CTL mapping
}
\usage{
CTLpermute(genotypes, phenotypes, geno.enc=c(1,2), pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, directory="permutations", saveFiles = FALSE, verbose=FALSE, \dots)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes (individuals x phenotypes) }
  \item{geno.enc}{ Specifying the encoding of the genotype matrix. }
  \item{pheno.col}{ Which phenotype column should we use in permutation analysis }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{n.perm}{ Number of permutations }
  \item{n.cores}{ Number of CPU cores to use during the analysis }
  \item{directory}{ Directory in which to save the output }
  \item{saveFiles}{ Save the permjutation matrices to HDD }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link{CTLscan}} function }
}
\details{
Permutation is one of the fundamental approach to determine probability when nothing is known about the underlying NULL-ditribution
\itemize{
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system
}
}
\value{
CTLpermute object, a list holding CTLscan objects ontained at each permutation. To view significance thresholds use 
the \code{\link{print.CTLpermute}} function.
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
  ctl_perms <- CTLpermute(ath.metab$genotypes, ath.metab$phenotypes, pheno.col = 1)
  ctl_perms
}

\keyword{methods}
