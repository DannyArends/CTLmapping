\name{plotTraits}
\alias{plotTraits}
\title{plotTraits - Trait vs Trait scatterplot, colored by the selected genetic locus }
\description{
  Trait vs Trait scatterplot, colored by the selected genetic locus
}
\usage{
plotTraits(genotypes, phenotypes, phenocol = c(1, 2), marker = 1, doRank = FALSE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{phenocol}{    Which phenotype column(s) should be plotted against each other, Default: phenotype 1 versus 2 }
  \item{marker}{      Which marker (column in genotypes) should be used to add genotype as a color of the dots. }
  \item{doRank}{      Transform quantitative data into ranked data before analyzing the slope. }
}
\details{
  TODO
}
\value{
  TODO
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
  data(ath.metabolites)                 # Arabidopsis Thaliana data set

  plotTraits(ath.metab$genotypes, ath.metab$phenotypes, marker=75, doRank = TRUE)
}
\keyword{methods}
