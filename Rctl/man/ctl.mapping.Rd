\name{CTLmapping}
\alias{CTLmapping}
\title{CTLmapping - Scan for correlated trait loci (CTL) }
\description{
  Scan for correlated trait loci (CTL)
}
\usage{
CTLmapping(genotypes, phenotypes, geno.enc=c(1,2), pheno.col = 1, method = c("pearson", "kendall", "spearman"), verbose = FALSE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes. (individuals x phenotypes) }
  \item{geno.enc}{ Specifying the encoding of the genotype matrix. }
  \item{pheno.col}{ Which phenotype column should we analyse. NOTE: CTLmapping is limited to a single trait, for multiple traits use \code{\link{CTLscan}}). }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{verbose}{ Be verbose. }
}
\details{
TODO
\itemize{
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system
}
}
\value{
CTLscan object, a list with at each index a CTL matrix (Rows: Phenotypes, 
Columns: Genetic markers) for the phenotype.
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
  ctl_result <- CTLmapping(ath.metab$genotypes, ath.metab$phenotypes)
  ctl_result
  plot(ctl_result)   #Plot a phenotype
}

\keyword{methods}
