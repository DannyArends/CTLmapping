\name{QCLscan.network}
\alias{QCLscan.network}

\title{QCLscan.network - Scan for quantitative correlation loci using a start phenotype }
\description{
	Scan for quantitative correlation loci using a start phenotype
}
\usage{
QCLscan.network(genotypes, phenotypes, pheno.col, qcl.threshold = 0.4, max.depth=10, verbose = TRUE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column should we analyse }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant. (See \code{\link{QCLpermute}}) }
  \item{max.depth}{ Number of times to adding new phenotypes which show QCL }
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
  qcl_result <- QCLscan.network(ath.metab$genotypes, ath.metab$phenotypes,pheno.col=2)
  qcl_result
}

\keyword{methods}
