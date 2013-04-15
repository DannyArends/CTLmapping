\name{scanSD}
\alias{scanSD}
\title{scanSD - Analyze the differences in Standard Deviation between genotypes between two traits }
\description{
  Analyze the differences in Standard Deviation between genotypes between two traits
}
\usage{
scanSD(genotypes, phenotypes, pheno.col=c(1,2), doRank = FALSE)
}

\arguments{
  \item{genotypes}{   Matrix of genotypes. (individuals x markers) }
  \item{phenotypes}{  Matrix of phenotypes. (individuals x phenotypes) }
  \item{pheno.col}{   Which phenotype column(s) should we analyse. Default analyse all phenotypes. }
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
  data(multitrait)      # Arabidopsis Thaliana (R/qtl cross object)

  sds <- scanSD(pull.geno(multitrait),pull.pheno(multitrait))
}

\keyword{methods}
