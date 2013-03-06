\name{plot.CTLscan}
\alias{plot.CTLscan}

\title{Plot CTL curves}

\description{
  Plot the CTL curve for a genome scan (the output of \code{\link{CTLscan}}).
}

\usage{
\method{plot}{CTLscan}(x, mapinfo = NULL, type = c("barplot","gwas","line"), onlySignificant = TRUE, significance = 0.05, plot.cutoff = FALSE, do.legend=TRUE, cex.legend=1.0, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{ mapinfo object, a matrix with in each row a marker information line with at least: Chromosome information in column 1, and cMorgan locations in column 2. }
  \item{type}{ Type of plot: Summed barplot, GWAS style plot or a basic line plot. }
  \item{onlySignificant}{ Plot only the significant contributions to the CTL profile. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{plot.cutoff}{ Adds a line at -log10(significance) and adds a legend showing the significance level. }
  \item{do.legend}{ Adds a legend showing which phenotypes contribute to the CTL profile. }
  \item{cex.legend}{ Maginification of the text in the legend. }
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
  data(ath.result)       # Arabidopsis thaliana results
  data(ath.metabolites)  # Arabidopsis thaliana data (phenotypes, genotypes and mapinfo

  plot(ath.result[[3]])
  plot(ath.result[[2]], mapinfo = ath.metab[[3]])
  plot(ath.result[[1]], mapinfo = ath.metab[[3]])
  plot(ath.result[[3]], mapinfo = ath.metab[[3]])
  plot(ath.result[[3]], mapinfo = ath.metab[[3]], type="gwas")
  plot(ath.result[[3]], mapinfo = ath.metab[[3]], type="line")
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
