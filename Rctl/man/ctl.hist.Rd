\name{hist.CTLobject}
\alias{hist.CTLobject}

\title{Plot histogram of CTL permutations}

\description{
  Plot histogram of CTL permutations (the output of \code{\link{CTLscan}}).
}

\usage{
\method{hist}{CTLobject}(x, pheno.col=1, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{pheno.col}{ Which phenotype column(s) should we analyse. Defaults to analyse all phenotype columns }
  \item{\dots}{ Passed to the function \code{\link[graphics]{image}} when it is called. }
}

\details{
  None.
}

\value{
  For a detailed description, see \code{\link{CTLprofiles}}
}

\examples{
  library(ctl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  ctl_result <- CTLscan(ath.metab$genotypes, ath.metab$phenotypes, 1:3, n.perm=50)
  hist(ctl_result,1:3)
}

\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Scan for CTL
    \item \code{\link{CTLprofiles}} - Extract CTL interaction profiles
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