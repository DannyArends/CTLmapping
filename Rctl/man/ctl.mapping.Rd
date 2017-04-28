\name{CTLmapping}
\alias{CTLmapping}
\title{CTLmapping - Scan for correlated trait loci (CTL) }
\description{
  Scan for correlated trait loci (CTL)
}
\usage{
CTLmapping(genotypes, phenotypes, phenocol = 1, nperm = 100, nthreads = 1,
strategy = c("Exact", "Full", "Pairwise"), adjust = TRUE, qtl = TRUE, verbose = FALSE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{phenocol}{    Which phenotype column(s) should we analyse. Default: Analyse a single phenotype. }
  \item{nperm}{       Number of permutations to perform. This parameter is not used when method="Exact". }
  \item{nthreads}{    Number of CPU cores to use during the analysis. }
  \item{strategy}{    The permutation strategy to use, either 
                      \itemize{
                        \item Exact: Uses exact calculations to calculate the likelihood of a difference 
                              in correlation: Cor(AA) - Cor(BB). Using a Bonferroni correction.
                        \item Full: Most powerful analysis method - Compensate for marker and trait 
                              correlation structure (Breitling et al.).
                        \item Pairwise: Suitable when we have a lot of markers and only a few 
                              traits (< 50) (human GWAS)- Compensates only for marker correlation structure.
                      }
                      Note: Exact is the default and fastest option it uses a normal distribution for 
                      estimating p-values and uses bonferoni correction. It has however the least power to 
                      detect CTLs, the two other methods (Full and Pairwise) perform permutations to
                      assign significance. }
  \item{adjust}{      Adjust p-values for multiple testing (only used when strategy = Exact). }
  \item{qtl}{         Use the internal slow QTL mapping method to map QTLs. }
  \item{verbose}{     Be verbose. }
}
\details{
TODO
\itemize{
  \item NOTE: Main bottleneck of the algorithm is the RAM available to the system
}
}
\value{
CTLscan, a list of:
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
\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Main function to scan for CTL
    \item \code{\link{CTLscan.cross}} - Use an R/qtl cross object with \code{\link{CTLscan}}
    \item \code{\link{CTLsignificant}} - Significant interactions from a \code{\link{CTLscan}}
    \item \code{\link{plot.CTLscan}} - Plot the CTL curve for a single trait
  }
}
\examples{
  library(ctl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  singlescan <- CTLmapping(ath.metab$genotypes, ath.metab$phenotypes, phenocol = 23)

  plot(singlescan)      # Plot the results of the CTL scan for the phenotype

  summary <- CTLsignificant(singlescan)
  summary               # Get a list of significant CTLs
}
\keyword{methods}
