\name{CTLscan.cross}
\alias{CTLscan.cross}
\title{CTLscan.cross - Scan for Correlated Trait Locus (CTL) (R/qtl cross object) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations (using an R/qtl cross object)
}
\usage{
CTLscan.cross(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=100, n.cores=2, genotype.values=c(1,2), directory="permutations", saveFiles = FALSE, verbose = FALSE)
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
  if(has_rqtl()){
    require(qtl)
    data(multitrait) # Arabidopsis Thaliana (R/qtl cross object)
    ctl_result <- CTLscan.cross(multitrait)
    ctl_result
    plot(ctl_result[[16]])   #Plot a phenotype
  }
}

\keyword{methods}