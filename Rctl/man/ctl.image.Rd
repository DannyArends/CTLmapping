\name{image.CTLobject}
\alias{image.CTLobject}
\title{Plot genome-wide CTL on multiple traits}
\description{
  Plot the CTL for genome-wide CTL on multiple traits (the output of \code{\link{CTLscan}}).
}
\usage{
\method{image}{CTLobject}(x, against = c("markers","phenotypes"), significance = 0.05,
col=whiteblack, do.grid=TRUE, grid.col = "white", verbose = FALSE, add=FALSE, 
breaks = c(0, 1, 2, 3, 10, 10000), marker_info, \dots)
}
\arguments{
  \item{x}{             An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{against}{       Plot which interaction matrice, options are: markers: the phenotype*marker or 
                        phenotypes: the phenotype*phenotypes matrix. }
  \item{significance}{  Significance threshold to set a genome wide False Discovery Rate (FDR). }
  \item{col}{           Color-range used in plotting. }
  \item{do.grid}{       When TRUE, grid lines are added to the plot. }
  \item{grid.col}{      Color used for the grid lines, only used when do.grid = TRUE. }
  \item{verbose}{       Be verbose. }
  \item{add}{           Add this plot to a previously opened plot window. }
  \item{breaks}{        See par. }
  \item{marker_info}{   Information used to plot chromosome lines. }
  \item{\dots}{         Passed to the function \code{\link[graphics]{image}} when it is called. }
}
\details{
  None.
}
\value{
  For a detailed description, see \code{\link{CTLprofiles}}
}
\examples{
  library(ctl)
  data(ath.result)       # Arabidopsis Thaliana results
  #Phenotype to phenotype matrix
  p2p_matrix <- image(ath.result, against="phenotypes")
  #Phenotype to marker matrix
  p2m_matrix <- image(ath.result, against="markers")
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
