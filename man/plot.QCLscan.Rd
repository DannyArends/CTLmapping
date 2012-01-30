\name{plot.QCLscan}
\alias{plot.QCLscan}

\title{Plot QCL curves}

\description{
  Plot the QCL curve for a genome scan (the output of \code{\link{QCLscan}}).
}

\usage{
\method{plot}{QCLscan}(x, addQTL = TRUE, onlySignificant = TRUE, qcl.threshold =0.6, do.legend=TRUE, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{addQTL}{ Plot the QTL on top of the QTL plot. }
  \item{onlySignificant}{ Plot only the significant contributions to the QCL profile. }
  \item{qcl.threshold}{ Threshold at which an QCL is considered significant, this is 
  used when onlySignificant == TRUE and no permutations are present in the QCLscan object. (See \code{\link{QCLscan}}). }
  \item{do.legend}{ Adds a legend showing which phenotypes contribute to the QCL profile. }
  \item{\dots}{ Passed to the function \code{\link[graphics]{plot}} when it is called. }
}
\value{
  None.
}

\details{
  None.
}

\examples{
  library(qcl)
  data(ath.result)       # Arabidopsis Thaliana dataset
  plot(ath.result[[1]])  # First phenotype
}

\seealso{
  \itemize{
    \item \code{\link{QCLscan}} - Scan for QCL
    \item \code{\link{print.QCLscan}} - Print a summary of a QCLscan
    \item \code{\link[graphics]{par}} - Plot parameters
    \item \code{\link[grDevices]{colors}} - Colors used in plotting
  }
}

\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\keyword{hplot}
