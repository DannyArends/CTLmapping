\name{ctl.lineplot}
\alias{ctl.lineplot}

\title{Lineplot CTL on multiple traits}

\description{
  Plot the CTL for genome-wide CTL on multiple traits (the output of \code{\link{CTLscan}}).
}

\usage{
ctl.lineplot(CTLobject, mapinfo, pheno.col, significance = 0.05, gap = 50, col = "orange", bg.col = "lightgray", cex = 1, verbose = FALSE)
}

\arguments{
  \item{CTLobject}{     An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{mapinfo}{       The mapinfo matrix with 3 columns: "Chr" - the chromosome number, "cM" - the 
                        location of the marker in centiMorgans and the 3rd column "Mbp" - The location 
                        of the marker in Mega basepairs. If supplied the marker names (rownames) 
                        should match those in the CTLobject (only significant markers will be annotated). }
  \item{pheno.col}{     Which phenotype results to plot. Defaults to plot all phenotypes. } 
  \item{significance}{  Significance threshold to set a genome wide False Discovery Rate (FDR). }
  \item{gap}{           The gap between chromosomes in cM. }
  \item{col}{           Line color used. }
  \item{bg.col}{        Node background color. }
  \item{cex}{           Global magnificantion factor for the image elements. }
  \item{verbose}{       Be verbose. }
}

\details{
  None.
}

\value{
  None.
}

\examples{
  require(ctl)
  data(ath.result)       # Arabidopsis Thaliana results
  data(ath.metabolites)  # Arabidopsis Thaliana data set

  todo <- c(1,3,4,5,6,8,9,10,11,12,14,17,18,19,22,23)
  op   <- par(mfrow = c(4,4))
  op   <- par(oma = c(0.1,0.1,0.1,0.1))
  op   <- par(mai = c(0.1,0.1,0.1,0.1))
  for(x in todo){ # Overview of the 16 traits with CTLs
    ctl.lineplot(ath.result, ath.metab$map, x, sign=0.1)
  }
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
