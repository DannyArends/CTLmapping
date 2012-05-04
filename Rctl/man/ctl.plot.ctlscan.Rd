\name{plot.CTLscan}
\alias{plot.CTLscan}

\title{Plot CTL curves}

\description{
  Plot the CTL curve for a genome scan (the output of \code{\link{CTLscan}}).
}

\usage{
\method{plot}{CTLscan}(x, addQTL = TRUE, onlySignificant = TRUE, significance= 0.05, do.legend = TRUE, cex.legend = 1, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{addQTL}{ Plot the QTL on top of the QTL plot. }
  \item{onlySignificant}{ Plot only the significant contributions to the CTL profile. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{do.legend}{ Adds a legend showing which phenotypes contribute to the CTL profile. }
  \item{cex.legend}{ Size of the text in the legend. }
  \item{\dots}{ Passed to the function \code{\link[graphics]{plot}} when it is called. }
}
\value{
  None.
}

\details{
  None.
}

\examples{
  library(ctl)
  data(ath.result)       # Arabidopsis Thaliana dataset
  plot(ath.result[[1]])  # First phenotype
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
