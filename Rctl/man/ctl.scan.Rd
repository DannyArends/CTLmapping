\name{CTLscan}
\alias{CTLscan}
\title{CTLscan - Scan for Correlated Trait Locus (CTL) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations
}
\usage{
CTLscan(genotypes, phenotypes, geno.enc=c(1,2), pheno.col = 1:ncol(phenotypes), method = c("pearson", "pearsonordered", "kendall", "spearman"), conditions = NULL, 
        n.perm=100, n.cores=2, directory="permutations", saveFiles = FALSE, verbose = FALSE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes. (individuals x phenotypes) }
  \item{geno.enc}{ Specifying the encoding of the genotype matrix. }
  \item{pheno.col}{ Which phenotype column(s) should we analyse. Defaults to analyse all phenotype columns. }
  \item{method}{ Character string indicating which correlation coefficient is to be computed. One of "pearson" (default), "kendall", or "spearman". }
  \item{conditions}{ A vector of experimental conditions applied during the experiment. These conditions will be used as covariates in the QTL modeling step
                     (NOT used when parameter have.qtl is supplied). }
  \item{n.perm}{ Number of permutations to perform. }
  \item{n.cores}{ Number of CPU cores to use during the analysis. }
  \item{directory}{ Directory in which to save the output. }
  \item{saveFiles}{ Save the permutation matrices to HDD. }
  \item{verbose}{ Be verbose. }
}
\details{
TODO
\itemize{
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system
}
}
\value{
CTLobject, a list with at each index a CTLscan object:
\itemize{
  \item $dcor  - Matrix of differential correlation scores for each trait at each marker
  \item $perms - Vector of maximums per marker obtained during permutations
  \item $ctls  - Matrix of LOD scores for CTL likelihood
}
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
  data(ath.metabolites)                    # Arabidopsis Thaliana dataset
  ctl_result <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, pheno.col=1:5, n.perm=30)
  ctl_result
  plot(ctl_result,3)                       # Plot the third phenotype
  image(ctl_result, against="markers")     # PxM matrix
  image(ctl_result, against="phenotypes")  # PxP matrix
}

\keyword{methods}
