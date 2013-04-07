\name{plot.CTLscan}
\alias{plot.CTLscan}

\title{Plot CTL results as bar, line or GWAS plot. }

\description{
  Plot the CTL results for a genome scan (the output of \code{\link{CTLscan}}) as a barplot, curved line 
  or GWAS plot.
}

\usage{
\method{plot}{CTLscan}(x, mapinfo = NULL, type = c("barplot","gwas","line"), addQTL = TRUE, onlySignificant = TRUE, significance = 0.05, plot.cutoff = FALSE, do.legend=TRUE, cex.legend=1.0, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{ The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the 
    location of the marker in centiMorgans and the 3rd column "Mbp" - The location of the marker 
    in Mega basepairs. If supplied the marker names (rownames) should match those in the 
    CTLobject (only significant markers will be annotated). }
  \item{type}{ Specify the type of plot to create. }
  \item{addQTL}{ Plot the QTL on top of the QTL plot. }
  \item{onlySignificant}{ Plot only the significant contributions to the CTL profile. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{plot.cutoff}{ Adds lines to show genome wide significance thresholds. }
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
  CTLscan3 <- ath.result[[3]]
  plot(CTLscan3)         # Third phenotype
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
