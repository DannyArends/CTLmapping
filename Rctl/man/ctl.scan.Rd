\name{CTLscan}
\alias{CTLscan}
\title{CTLscan - Scan for Correlated Trait Locus (CTL) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations
}
\usage{
CTLscan(genotypes, phenotypes, pheno.col, n.perm=100, strategy = c("Exact", "Full", "Pairwise"), 
        conditions = NULL, qtls = NULL, n.cores = 1, geno.enc = c(1,2), verbose = FALSE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{pheno.col}{   Which phenotype column(s) should we analyse. Default analyse all phenotypes. }
  \item{n.perm}{      Number of permutations to perform. This parameter is not used when method="Exact". }
  \item{strategy}{    The permutation strategy to use, either 
                      \itemize{
                        \item Exact - Uses exact calculations to calculate the likelihood of a difference 
                              in correlation: Cor(AA) - Cor(BB). Using a Bonferroni correction.
                        \item Full - Most powerful analysis method - Compensate for marker and trait 
                              correlation structure (Breitling et al.).
                        \item Pairwise - Suitable when we have a lot of markers and only a few 
                              traits (< 50) (human GWAS)- Compensates only for marker correlation structure.
                      }
                      Note: Exact is the default and fastest option it uses a normal distribution for 
                      estimating p-values and uses bonferoni correction. It has however the least power to 
                      detect CTLs, the two other methods (Full and Pairwise) perform permutations to
                      assign significance. }
  \item{conditions}{  A vector of experimental conditions applied during the experiment. These conditions 
                      are used as covariates in the QTL modeling step. }
  \item{qtls}{        Used to provide QTL results, when external QTL results are available. }
  \item{n.cores}{     Number of CPU cores to use during the analysis. }
  \item{geno.enc}{    Specifying the encoding of the genotype matrix. }
  \item{verbose}{     Be verbose. }
}

\details{
By default the algorithm will not do QTL mapping, the qtl component of the output is an vector of 0 scores 
for LOD. This is to remove some computational burden, please use the have.qtls parameter to provide QTL data. 
Some computational bottleneck of the algorithm are:
\itemize{
  \item RAM available to the system with large number of markers (100K+) and/or  phenotypes (100K+).
  \item Computational time with large sample sizes (5000+) and/or huge amount of phenotype data (100K+).
  \item Very very huge amounts of genotype markers (1M+)
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

  ctlscan <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, pheno.col=1:4)
  ctlscan

  summary <- CTLsignificant(ctlscan)    # Matrix of Trait, Marker, Trait interactions 
  summary                               # Get a list of significant CTLs

  nodes <- ctl.lineplot(ctlscan, ath.metab$map)  # Line plot the phenotypes
  nodes
}

\keyword{methods}
