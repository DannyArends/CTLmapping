\name{QCLmapping}
\alias{QCLmapping}
\title{QCLmapping - Scan for quantitative correlation loci }
\description{
  Scan for quantitative correlation loci
}
\usage{
QCLmapping(genotypes, phenotypes, pheno.col = 1, method = c("pearson", "kendall", "spearman"), genotype.values=c(1,2), verbose = FALSE)
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
QCLscan object, a list with at each index a QCL matrix (Rows: Phenotypes, 
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
  library(qcl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  qcl_result <- QCLmapping(ath.metab$genotypes, ath.metab$phenotypes)
  qcl_result
  plot(qcl_result)   #Plot a phenotype
}

\keyword{methods}
