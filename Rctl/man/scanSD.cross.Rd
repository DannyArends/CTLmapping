\name{scanSD.cross}
\alias{scanSD.cross}
\title{scanSD.cross - Analyze the differences in standard deviation between two traits at a certain genetic marker (R/qtl cross object) }
\description{
  Analyze the differences in standard deviation between two traits at a certain genetic marker
}
\usage{
scanSD.cross(cross, pheno.col = c(1,2), doRank = FALSE)
}
\arguments{
  \item{cross}{       An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
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

  sds <- scanSD.cross(multitrait)
}
\keyword{methods}
