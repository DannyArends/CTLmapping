\name{QCLscan.cross}
\alias{QCLscan.cross}
\title{QCLscan.cross - Scan for quantitative correlation loci (R/qtl cross object) }
\description{
  Scan for Quantitative Correlation Loci (QCL) in inbred populations (using an R/qtl cross object)
}
\usage{
QCLscan.cross(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose = FALSE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype column should we analyse }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{n.perm}{ Number of permutations }
  \item{n.cores}{ Number of CPU cores to use during the analysis }
  \item{genotype.values}{ Two values specifying the encoding of the genotype matrix }
  \item{directory}{ Directory in which to save the output }
  \item{saveFiles}{ Save the permjutation matrices to HDD }
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
    plot(qcl_result[[16]])   #Plot a phenotype
  }
}

\keyword{methods}
