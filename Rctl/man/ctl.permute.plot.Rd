\name{plot.CTLpermute}
\alias{plot.CTLpermute}
\title{Differential correlation versus likelihood plotted in curves}
\description{
  Differential correlation versus likelihood plot curves.
}
\usage{
\method{plot}{CTLpermute}(x, type="s", \dots)
}
\arguments{
  \item{x}{       An object of class \code{"CTLscan"}. }
  \item{type}{    What type of plot should be drawn. for possible options see \code{\link[base]{plot}}. }
  \item{\dots}{   Passed to the function \code{\link[base]{plot}} when it is called. }
}
\value{
  None.
}
\details{
  None.
}
\examples{
  library(ctl)
  data(ath.result) # Arabidopsis Thaliana dataset
  plot(ath.result[[1]]$perms)
}
\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Scan for CTL
    \item \code{\link{print.CTLscan}} - Print a summary of a CTLscan
    \item \code{\link[graphics]{par}} - Plot parameters
    \item \code{\link[grDevices]{colors}} - Colors used in plotting
  }
}
\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\keyword{hplot}
