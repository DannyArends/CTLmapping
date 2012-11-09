\name{CTLmapping}
\alias{CTLmapping}
\title{CTLmapping - Scan for correlated trait loci (CTL) }
\description{
  Scan for correlated trait loci (CTL)
}
\usage{
CTLmapping(genotypes, phenotypes, pheno.col=1, method = c("Exact","Power","Adjacency"), n.perms=100, strategy = c("Full", "Pairwise"), have.qtls=NULL, geno.enc=c(1,2), verbose = FALSE)
}

\arguments{
  \item{genotypes}{ Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{ Matrix of phenotypes. (individuals x phenotypes) }
  \item{pheno.col}{ Which phenotype column should we analyse. NOTE: CTLmapping is limited to a single trait, for multiple traits use \code{\link{CTLscan}}). }
  \item{method}{ 
    We provide 3 ways of mapping correlation differences across the genome:
        Exact: Use Cor->Z transformation to calculate the likelihood of a difference in correlation: Cor(AA) - Cor(BB) 
        Power: Power analysis using the squared difference in correlation (Cor(AA) - Cor(BB))^2
        Adjacency: Power analysis using the squared difference in squared corrrelation (sign*Cor(AA)^2 - sign*Cor(BB)^2)^2
    Note: Exact is the default and fastest option it uses a normal distribution for estimating p-values and uses bonferoni correction, 
    the two other methods (Power and Adjacency) perform permutations to assign significance.
  }  
  \item{n.perms}{ Number of permutations to perform. }
  \item{strategy}{ The permutation strategy to use, either Full (Compensate for marker and trait correlation structure) or Pairwise (Compensate for 
    marker correlation structure). This parameter is not used when method="Exact". 
  }
  \item{have.qtls}{ Can be used to provide QTL results, when external QTL results are available. }
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
  singlescan <- CTLmapping(ath.metab$genotypes, ath.metab$phenotypes, pheno.col=23)

  plot(singlescan)      # Plot the results of the CTL scan for the phenotype

  summary <- CTLsignificant(singlescan)
  summary               # Get a list of significant CTLs
}

\keyword{methods}
