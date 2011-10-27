\name{QCLscan}
\alias{QCLscan}
\alias{scanQCL}
\title{QCLscan - Scan for quantitative correlation loci }
\description{
	Scan for quantitative correlation loci
}
\usage{
QCLscan(cross, pheno.col, verbose = FALSE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype column should we analyse }
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
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  library(qcl)
  data(multitrait)
  multitrait <- fill.geno(multitrait)
  qcl_result <- QCLscan(multitrait)
  qcl_result
  names(qcl_result)     #Names of the analysed phenotypes
  plot(qcl_result,16)   #Plot a phenotype
}

\keyword{methods}
