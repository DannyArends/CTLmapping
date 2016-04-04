\name{scanSlopes.cross}
\alias{scanSlopes.cross}
\title{scanSlopes.cross - Create a slope difference profile between two traits (R/qtl cross object) }
\description{
  Create a slope difference profile between two traits (using an R/qtl cross object)
}
\usage{
scanSlopes.cross(cross, phenocol = 1, doRank = FALSE, verbose = FALSE)
}

\arguments{
  \item{cross}{       An object of class \code{cross}. See \code{\link[qtl]{read.cross}} for details. }
  \item{phenocol}{    Which phenotype column(s) should we analyse. Default: Analyse phenotype 1 }
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
  data(multitrait)      # Arabidopsis Thaliana (R/qtl cross object)

  slopes <- scanSlopes.cross(multitrait)
  image(1:nrow(slopes), 1:ncol(slopes), -log10(slopes))
}

\keyword{methods}
