\name{CTLscan}
\alias{CTLscan}
\title{CTLscan - Scan for Correlated Trait Locus (CTL) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations
}
\usage{
CTLscan(genotypes, phenotypes, pheno.col, method = c("Exact","Power","Adjacency"), n.perms=100, strategy = c("Full", "Pairwise"), 
        conditions = NULL, have.qtls=NULL, n.cores = 1, geno.enc = c(1,2), verbose = FALSE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{pheno.col}{   Which phenotype column(s) should we analyse. Default analyse all phenotypes. }
  \item{method}{      We provide 3 ways of mapping correlation differences across the genome:
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
                      estimating p-values and uses bonferoni correction. It has however the least power to 
                      detect CTLs, the two other methods (Power and Adjacency) perform permutations to
                      assign significance. }
  \item{n.perms}{     Number of permutations to perform. This parameter is not used when method="Exact". }
  \item{strategy}{    The permutation strategy to use, either Full (Compensate for marker and trait 
                      correlation structure) or Pairwise (Compensate only for marker correlation
                      structure). This parameter is not used when method="Exact". }
  \item{conditions}{  A vector of experimental conditions applied during the experiment. These conditions 
                      are used as covariates in the QTL modeling step. }
  \item{have.qtls}{   Used to provide QTL results, when external QTL results are available. }
  \item{n.cores}{     Number of CPU cores to use during the analysis. }
  \item{geno.enc}{    Specifying the encoding of the genotype matrix. }
  \item{verbose}{     Be verbose. }
}

\details{
Some computational bottleneck of the algorithm are:
\itemize{
  \item RAM available to the system
  \item Computational time with large sample sizes and/or huge amount of phenotype data
  \item Huge amounts of genotype markers
}
Some way of avoiding these problems are: CTL mapping using only a single chromosome at a time and/or selecting a smaller subsets of phenotype data for analysis.
}

\value{
CTLobject, a list with at each index (i) an CTLscan object:
\itemize{
  \item $dcor  - Matrix of Z scores (method=Exact), or Power/Adjacency Z scores or for each trait at each
                 marker (n.markers x n.phenotypes)
  \item $perms - Vector of maximum scores obtained during permutations  (n.perms)
  \item $ctl   - Matrix of LOD scores for CTL likelihood of phenotype i (n.markers x n.phenotypes)
  \item $qtl   - Vector of LOD scores for QTL likelihood of phenotype i (n.markers)
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
    \item \code{\link{CTLscan.cross}} - Use an R/qtl cross object with \code{\link{CTLscan}}
    \item \code{\link{CTLsignificant}} - Significant interactions from a \code{\link{CTLscan}}
    \item \code{\link{CTLnetwork}} - Create a CTL network from a \code{\link{CTLscan}}
    \item \code{\link{image.CTLobject}} - Heatmap overview of a \code{\link{CTLscan}}
    \item \code{\link{plot.CTLscan}} - Plot the CTL curve for a single trait
  }
}
\examples{
  library(ctl)
  data(ath.metabolites)                 # Arabidopsis Thaliana data set

  ctlscan <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, pheno.col=1:4, method="Power")
  ctlscan

  summary <- CTLsignificant(ctlscan)    # Matrix of Trait, Marker, Trait interactions 
  summary                               # Get a list of significant CTLs

  nodes <- ctl.lineplot(ctlscan, ath.metab$map)  # Line plot the phenotypes
  nodes
}

\keyword{methods}
