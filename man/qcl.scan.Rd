\name{QCL.scan}
\alias{QCL.scan}
\title{QCL.scan - Scan for quantitative correlation loci }
\description{
	Scan for quantitative correlation loci
}
\usage{
QCL.scan(cross, pheno.col, verbose = FALSE)
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
QCLscan object, a list with at each index a QCL matrix (Rows: Phenotypes, Columns: Genetic markers) for the phenotype.
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
  require(qtl)
  data(multitrait)
  multitrait = fill.geno(multitrait)
  res <- QCL.scan(multitrait)
  dim(res)              #show the dimensions
  attr(res[[1]],"name") #show the name of the first phenotype
  res[[1]][1:10,1:10]   #A piece of the QCL matrix
}

\keyword{methods}
