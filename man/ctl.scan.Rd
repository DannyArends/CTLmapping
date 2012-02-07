\name{CTLscan}
\alias{CTLscan}
\title{CTLscan - Scan for Correlated Trait Locus (CTL) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations
}
\usage{
CTLscan(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), method = c("pearson", "kendall", "spearman"), n.perm=0, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose = FALSE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column(s) should we analyse. Defaults to analyse all phenotype columns }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{n.perm}{ Number of permutations }
  \item{n.cores}{ Number of CPU cores to use during the analysis }
  \item{genotype.values}{ Two values specifying the encoding of the genotype matrix }
  \item{directory}{ Directory in which to save the output }
  \item{saveFiles}{ Save the permjutation matrices to HDD }
  \item{verbose}{ Be verbose }
}
\details{
TODO
\itemize{
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system
}
}
\value{
CTL object, a list with at each index a $ctlscan, $permutations, $lodmatrix
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
  ctl_result <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, 1:3, n.perm=20)
  ctl_result
  plot(ctl_result)   #Plot a phenotype
  image(ctl_result, against="markers")     #pxm matrix
  image(ctl_result, against="phenotypes")  #pxp matrix
}

\keyword{methods}
