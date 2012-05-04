\name{image.CTLobject}
\alias{image.CTLobject}

\title{Plot genome-wide CTL on multiple traits}

\description{
  Plot the CTL for genome-wide CTL on multiple traits (the output of \code{\link{CTLscan}}).
}

\usage{
\method{image}{CTLobject}(x, against = c("markers","phenotypes"), onlySignificant = FALSE, significance = 0.05, 
      col=c("white",gray.colors(10)[10:1]), do.grid=TRUE, grid.col = "black", verbose = FALSE, add=FALSE, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{against}{ Which interaction matrice to image, the phenotype*marker or the phenotype*phenotypes, 
  this parameter is passed to the \code{\link{CTLprofiles}} function. }
  \item{onlySignificant}{ Plot only the significant contributions to the CTL profile. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{col}{ Color-range used in plotting. }
  \item{do.grid}{ Draw grid lines. }
  \item{grid.col}{ Colors used to draw the grid. }
  \item{verbose}{ Be verbose. }
  \item{add}{ Add this plot to a previously opened plot window. }
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
  data(ath.result)       # Arabidopsis Thaliana dataset
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
