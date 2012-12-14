\name{ctl-package}
\alias{ctl}
\alias{CTL}
\docType{package}
\title{
  CTL - CTL mapping in experimental crosses
}
\description{
  Analysis of experimental crosses to identify genetic markers associated with correlation 
  changes in quantitative traits (CTL). The additional correlation information obtained 
  can be combined with QTL information to perform de novo reconstruction of interaction 
  networks.

  For more background information about the method we refer to the methodology article 
  published in XX (201X).

  The R package is a basic iomplementation and it includes the following core functionality:
  \itemize{
    \item \code{\link{CTLscan}} - Main function to scan for CTL.
    \item \code{\link{CTLsignificant}} - Significant interactions from a \code{\link{CTLscan}}.
    \item \code{\link{CTLnetwork}} - Create a CTL network from a \code{\link{CTLscan}}.
    \item \code{\link{image.CTLobject}} - Heatmap overview of a CTLscan.
    \item \code{\link{plot.CTLscan}} - Plot the CTL curve for a single trait.
    \item \code{\link{ctl.circle}} - Circle plot CTLs on single and multiple traits.
    \item \code{\link{ctl.lineplot}} - Line plot CTLs on single and multiple traits.
    \item \code{\link{CTLprofiles}} - Extract CTL interaction profiles.
  }
  For all these functions we also provide examples and demonstrations on real genetical 
  genomics data. We thank all contributors for publishing their data online and will 
  accept submissions of intrestion datasets, currently ctl provides:
  \itemize{
    \item \code{\link{ath.metabolites}} - Metabolite expression data from Arabidopsis Thaliana
    \item \code{\link{ath.churchill}} - Metabolite expression data from Arabidopsis Thaliana
    \item \code{\link{yeast.brem}} - Gene expression data from Saccharomyces cerevisiae
  }
}
\details{
  More detailed information and/or examples are given per function as needed.
  Some additional functionality:
  \itemize{
    \item \code{\link{basic.qc}} - Some basic quality checks for phenotype and genotype data
    \item \code{\link{CTLscan.cross}} - Use an R/qtl cross object with \code{\link{CTLscan}}
  }
}
\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com}\cr
  Contributions from: Bruno Tesson, Pjotr Prins and Ritsert C. Jansen
}
\references{
  \itemize{
    \item TODO
  }
}
\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Scan for CTL
    \item \code{\link{CTLscan.cross}} - Use an R/qtl cross object with \code{\link{CTLscan}}
  }
}
\keyword{package}
