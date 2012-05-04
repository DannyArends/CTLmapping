\name{CTLpermute.cross}
\alias{CTLpermute.cross}
\title{CTLpermute.cross - Permutation routine for CTL mapping (R/qtl cross object) }
\description{
  Permutation routine for CTL mapping (using an R/qtl cross object)
}
\usage{
CTLpermute.cross(cross, pheno.col, method = c("pearson", "kendall", "spearman"), n.perm=10, n.cores=2, directory="permutations", saveFiles=FALSE, verbose=FALSE, \dots)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{pheno.col}{ Which phenotype column should we use in permutation analysis }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{n.perm}{ Number of permutations }
  \item{n.cores}{ Number of CPU cores to use during the analysis }
  \item{directory}{ Directory in which to save the output }
  \item{saveFiles}{ Save the permjutation matrices to HDD }
  \item{verbose}{ Be verbose }
  \item{\dots}{ Passed to the underlying \code{\link{CTLscan}} function }
}
\details{
TODO
\itemize{
  \item NOTE: In general using the interface function \code{\link{CTLscan}} is advised.
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system.
}
}
\value{
CTLpermute object, a list holding CTLscan objects ontained at each permutation. To view significance 
thresholds at different significance levels use the \code{\link{print.CTLpermute}} function.
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
  \dontrun{
    require(qtl)
    data(multitrait) # Arabidopsis Thaliana (R/qtl cross object)
    ctl_perms <- CTLpermute.cross(multitrait,1)
    ctl_perms
  }
}

\keyword{methods}
