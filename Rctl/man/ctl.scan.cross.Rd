\name{CTLscan.cross}
\alias{CTLscan.cross}
\title{CTLscan.cross - Scan for Correlated Trait Locus (CTL) (R/qtl cross object) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations (using an R/qtl cross object)
}
\usage{
CTLscan.cross(cross, ...)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{...}{ Passed to \code{\link{CTLscan}} function: 
    \itemize{
      \item phenocol   - Which phenotype column should we analyse.
      \item method     - We provide 3 ways of mapping correlation differences across the genome:
                         \itemize{
                           \item Exact: Uses a Correlation to Z score transformation to calculate the
                                 likelihood of a difference in correlation: 
                                 Cor(AA) - Cor(BB)
                           \item Power: More powerful analysis method using the squared difference in     
                                 correlation:
                                 (Cor(AA) - Cor(BB))^2
                           \item Adjacency: Adjacency method which using the squared difference in squared   
                                 correlation, but keeping the sign of correlation:
                                 (sign*Cor(AA)^2 - sign*Cor(BB)^2)^2
                         }
                         Note: Exact is the default and fastest option it uses a normal distribution for 
                         estimating p-values and uses bonferoni correction. It has however the least power 
                         to detect CTLs, the two other methods (Power and Adjacency) perform permutations to
                         assign significance.
      \item n.perm     - Number of permutations to perform.
      \item strategy   - The permutation strategy to use, either Full (Compensate for marker and trait
                         correlation structure) or Pairwise (Compensate for marker correlation structure). 
                         This parameter is not used when method="Exact".
      \item conditions - A vector of experimental conditions applied during the experiment. These conditions 
                         will be used as covariates in the QTL modeling step.
      \item qtls       - Can be used to provide QTL results, when external QTL results are available.
      \item n.cores    - Number of CPU cores to use during the analysis.
      \item verbose    - Be verbose.
    }
  }
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
\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Main function to scan for CTL
    \item \code{\link{CTLsignificant}} - Significant interactions from a \code{\link{CTLscan}}
    \item \code{\link{CTLnetwork}} - Create a CTL network from a \code{\link{CTLscan}}
    \item \code{\link{image.CTLobject}} - Heatmap overview of a \code{\link{CTLscan}}
    \item \code{\link{plot.CTLscan}} - Plot the CTL curve for a single trait
  }
}

\examples{
  library(ctl)
  data(multitrait)      # Arabidopsis Thaliana (R/qtl cross object)

  mtrait <- calc.genoprob(multitrait)          # Calculate genotype probabilities
  qtls   <- scanone(mtrait, pheno.col = 1:4)   # Scan for QTLS using R/qtl
  map    <- qtls[,1:2]                         # Get the genetic map
  qtls   <- qtls[,-(1:2)]                      # Get the LOD scores

  ctls   <- CTLscan.cross(mtrait, phenocol = 1:4, qtls = qtls)
  
  ctl.lineplot(ctls, map, significance = 0.05) # Line plot all the phenotypes

  summary <- CTLsignificant(ctls)              # Get a list of significant CTLs
  summary
}

\keyword{methods}
