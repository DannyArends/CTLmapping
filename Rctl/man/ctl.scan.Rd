\name{CTLscan}
\alias{CTLscan}
\title{CTLscan - Scan for Correlated Trait Locus (CTL) }
\description{
  Scan for Correlated Trait Locus (CTL) in populations
}
\usage{
CTLscan(genotypes, phenotypes, pheno.col = 1:ncol(phenotypes), method = c("Exact","Power","Adjacency"), n.perms=100, strategy = c("Full", "Pairwise"), 
        conditions = NULL, have.qtls=NULL, n.cores = 1, geno.enc = c(1,2), verbose = FALSE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes. (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column(s) should we analyse. Defaults to analyse all phenotype columns. }
  \item{method}{ 
    We provide 3 ways of mapping correlation differences across the genome:
        Exact: Use Cor->Z transformation to calculate the likelihood of a difference in correlation: Cor(AA) - Cor(BB) 
        Power: Power analysis using the squared difference in correlation (Cor(AA) - Cor(BB))^2
        Adjacency: Power analysis using the squared difference in squared corrrelation (sign*Cor(AA)^2 - sign*Cor(BB)^2)^2
    Note: Exact is the default and fastest option it uses a normal distribution for estimating p-values and uses bonferoni correction, 
    the two other methods (Power and Adjacency) perform permutations to assign significance.
  }
  \item{n.perms}{ Number of permutations to perform. This parameter is not used when method="Exact". }
  \item{strategy}{ The permutation strategy to use, either Full (Compensate for marker and trait correlation structure) or Pairwise (Compensate for 
    marker correlation structure). This parameter is not used when method="Exact". 
  }
  \item{conditions}{ A vector of experimental conditions applied during the experiment. These conditions are used as covariates in the QTL modeling 
    step. 
  }
  \item{have.qtls}{ Can be used to provide QTL results, when external QTL results are available. }
  \item{n.cores}{ Number of CPU cores to use during the analysis. }
  \item{geno.enc}{ Specifying the encoding of the genotype matrix. }
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
  data(ath.metabolites)                    # Arabidopsis Thaliana dataset

  ctl_result <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, pheno.col=1:3)
  ctl_result
  as.numeric(  ctl_result[[1]]$ctl[,1])
  as.numeric(  ctl_result[[1]]$dcor[,1])

  summary <- CTLsignificant(ctl_result)    # Matrix of Trait, Marker, Trait interactions 
  summary                                  # Get a list of significant CTLs

  plot(ctl_result,3)                       # Plot the third phenotype
  image(ctl_result, against="markers")     # PxM matrix
  image(ctl_result, against="phenotypes")  # PxP matrix
}

\keyword{methods}
