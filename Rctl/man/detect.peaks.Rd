\name{detect.peaks}
\alias{detect.peaks}
\title{detect.peaks - Peak detection algorithm to 'flatten' data above a certain threshold }
\description{
  Peak detection algorithm to 'flatten' data above a certain threshold
}
\usage{
detect.peaks(data, chrEdges = c(1), threshold = 4, verbose = FALSE)
}
\arguments{
  \item{data}{        A vector of scores per marker/locus. }
  \item{chrEdges}{    Start positions of the chromosomes. }
  \item{threshold}{   Threshold to determine regions. }
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
  #TODO
}
\keyword{methods}
