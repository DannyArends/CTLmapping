\name{QCLscan.network.cross}
\alias{QCLscan.network.cross}

\title{QCLscan.network.cross - Scan for quantitative correlation loci using a start phenotype }
\description{
	Scan for quantitative correlation loci using a start phenotype  (using an R/qtl cross object)
}
\usage{
QCLscan.network.cross(cross, pheno.col, qcl.threshold = 0.4, max.depth=10, verbose = TRUE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
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
  \dontrun{
    require(qtl)
    data(multitrait) # Arabidopsis Thaliana (R/qtl cross object)
    qcl_result <- QCLscan.network.cross(multitrait,pheno.col=2)
    qcl_result
  }
}

\keyword{methods}
