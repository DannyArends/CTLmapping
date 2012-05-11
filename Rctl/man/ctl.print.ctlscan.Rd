\name{print.CTLscan}
\alias{print.CTLscan}

\title{Print the results of a single phenotype CTL scan}

\description{
  Print the results of a single phenotype CTL scan produced by either \code{\link{CTLmapping}} (Single phenotype scan) or 
  \code{\link{CTLscan}} (Multi phenotype scan).
}

\usage{
\method{print}{CTLscan}(x, \dots)
}

\arguments{
 \item{x}{ An object of class \code{"CTLscan"}. This is a single element from an \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
 \item{\dots}{ Ignored.} 
}

\details{
  None.
}

\value{
  None. 
}

\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\examples{
  #TODO
}

\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Scan for CTL
    \item \code{\link{plot.CTLscan}} - Plot a CTLscan object
  }
}

\references{
  TODO
}

\keyword{print}
