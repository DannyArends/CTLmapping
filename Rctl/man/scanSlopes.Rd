\name{scanSlopes}
\alias{scanSlopes}
\title{scanSlopes - Create a slope difference profile between two traits }
\description{
  Create a slope difference profile between two traits
}
\usage{
scanSlopes(genotypes, phenotypes, phenocol = 1, doRank = FALSE, verbose = FALSE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{phenocol}{    Which phenotype column(s) should we analyse. Default: Analyse phenotype 1. }
  \item{doRank}{      Transform quantitative data into ranked data before analyzing the slope. }
  \item{verbose}{     Be verbose. }
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

  slopes <- scanSlopes(ath.metab$genotypes, ath.metab$phenotypes, phenocol = 2)
  image(1:nrow(slopes), 1:ncol(slopes), -log10(slopes))
}

\keyword{methods}
