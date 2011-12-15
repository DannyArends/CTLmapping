\name{QCLscan}
\alias{QCLscan}
\title{QCLscan - Scan for quantitative correlation loci }
\description{
	Scan for quantitative correlation loci
}
\usage{
QCLscan(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), method = c("pearson", "kendall", "spearman"), n.perm=0, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose = FALSE)
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
QCL object, a list with at each index a $qclscan, $permutations, $lodmatrix
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
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes, 1:3, n.perm=20)
  qcl_result
  plot(qcl_result)   #Plot a phenotype
  image(qcl_result, against="markers")     #pxm matrix
  image(qcl_result, against="phenotypes")  #pxp matrix
}

\keyword{methods}
