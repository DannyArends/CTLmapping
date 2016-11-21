\name{QTLmapping}
\alias{QTLmapping}
\title{QTLmapping - QTL mapping method for CTL analysis }
\description{
  Internal QTL mapping method used by the CTL analysis, associates every column in the genotypes with a single phenotype
}
\usage{
QTLmapping(genotypes, phenotypes, phenocol = 1, verbose = TRUE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{phenocol}{    Which phenotype column(s) should we analyse. Default: Analyse a single phenotype. }
  \item{verbose}{     Be verbose. }
}
\details{
TODO
\itemize{
  \item NOTE: Slow approach, it is adviced to use your own QTL mapping data
}
}
\value{
vector of LOD scores for each genotype column, for phenotype column phenocol
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
  qtldata <- QTLmapping(ath.metab$genotypes, ath.metab$phenotypes, phenocol = 23)
  plot(qtldata)         # Plot the results of the QTL scan for the phenotype
}
\keyword{methods}
