\name{QCLscan.cross}
\alias{QCLscan.cross}
\title{QCLscan.cross - Scan for quantitative correlation loci (R/qtl cross object) }
\description{
	Scan for quantitative correlation loci (using an R/qtl cross object)
}
\usage{
QCLscan.cross(cross, pheno.col, verbose = FALSE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype column should we analyse }
  \item{verbose}{ Be verbose. }
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
  if(has_rqtl()){
    require(qtl)
    data(multitrait) # Arabidopsis Thaliana (R/qtl cross object)
    qcl_result <- QCLscan.cross(multitrait)
    qcl_result
    plot(qcl_result,16)   #Plot a phenotype
  }
}

\keyword{methods}
