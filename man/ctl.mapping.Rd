\name{CTLmapping}
\alias{CTLmapping}
\title{CTLmapping - Scan for correlated trait loci (CTL) }
\description{
  Scan for correlated trait loci (CTL)
}
\usage{
CTLmapping(genotypes, phenotypes, pheno.col = 1, method = c("pearson", "kendall", "spearman"), genotype.values=c(1,2), verbose = FALSE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column should we analyse }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{genotype.values}{ Two values specifying the encoding of the genotype matrix }
  \item{verbose}{ Be verbose }
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
