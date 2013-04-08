\name{plot.CTLscan}
\alias{plot.CTLscan}
\title{Plot CTL results as bar, line or GWAS plot. }
\description{
  Plot the CTL results for a genome scan (the output of \code{\link{CTLscan}}) as a barplot, curved line 
  or GWAS plot.
}
\usage{
\method{plot}{CTLscan}(x, mapinfo, type = c("barplot","gwas","line"), onlySignificant = TRUE, 
                       significance = 0.05, gap = 25, plot.cutoff = FALSE, do.legend=TRUE, 
                       cex.legend=1.0, ydim=NULL, ylab="-log10(P-value)", \dots)
}
\arguments{
  \item{x}{                 An object of class \code{"CTLscan"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{           The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the 
                            location of the marker in centiMorgans and the 3rd column "Mbp" - The location 
                            of the marker in Mega basepairs. If supplied the marker names (rownames) 
                            should match those in the CTLobject. }
  \item{type}{              Type of plot: Summed barplot, GWAS style plot or a basic line plot. }
  \item{onlySignificant}{   Plot only the significant contributions to the CTL profile. }
  \item{significance}{      Significance threshold for setting a genomewide FDR. }
  \item{gap}{               Gap in Cm between chromosomes. }
  \item{plot.cutoff}{       Adds a line at -log10(significance) and adds a legend showing the 
                            significance level. }
  \item{do.legend}{         Adds a legend showing which phenotypes contribute to the CTL profile. }
  \item{cex.legend}{        Maginification of the text in the legend. }
  \item{ydim}{              Dimension of the y-axis, if NULL then it will be calculated. }
  \item{ylab}{              Label for the y-axis. }
  \item{\dots}{             Passed to the function \code{\link[graphics]{plot}} when it is called. }
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
