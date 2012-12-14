\name{ctl.circle}
\alias{ctl.circle}

\title{Circleplot CTL on multiple traits}

\description{
  Plot the CTL for genome-wide CTL on multiple traits (the output of \code{\link{CTLscan}}).
}

\usage{
ctl.circle(CTLobject, mapinfo, pheno.col, significance = 0.05, gap = 50, cex = 1, verbose = FALSE)
}

\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{ Genetic map chromosome and location information (Colums: Chr cMPos). }
  \item{pheno.col}{ Which phenotype results to plot. Defaults to plot all phenotypes. } 
  \item{significance}{ Significance threshold to set a genome wide False Discovery Rate (FDR). }
  \item{gap}{ Gap between chromosomes in cM. }
  \item{cex}{ Global magnificantion factor for the image elements. }
  \item{verbose}{ Be verbose. }
}

\details{
  None.
}

\value{
  None.
}

\examples{
  library(ctl)
  data(ath.result)       # Arabidopsis Thaliana results
  ctl.circle(ath.result, ath.metab$map, sign=0.001)
  ctl.circle(ath.result, ath.metab$map, pheno.col = 1:6, sign=0.01)
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
