\name{QCLscan}
\alias{QCLscan}
\title{QCLscan - Scan for quantitative correlation loci }
\description{
	Scan for quantitative correlation loci
}
\usage{
QCLscan(cross, marker.col, QCL.threshold=0.25, significant = 0, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, saveRdata=FALSE, snow=TRUE, verbose=TRUE)
}

\arguments{
  \item{cross}{ An object of class \code{cross}. See
                \code{\link[qtl]{read.cross}} for details. }
  \item{marker.col}{ Which marker should we analyse }
  \item{QCL.threshold}{ Changes in correlation^2 >= QCL.threshold }
  \item{significant}{ Used as a cut-off for plotting }
  \item{method}{ Method of calculating correlation }
  \item{directory}{ Directory in which to save the output }
  \item{doplot}{ Create per marker plots }
  \item{writefile}{ Write intermediate output files }
  \item{saveRdata}{ Save the correlation matrices as Rdata object to the directory }
  \item{snow}{ Use the SNOW library to use multiple cpus }
  \item{verbose}{ Be verbose }
}
\details{
TODO
\itemize{
  \item NOTE: QCLscan is optimized for multicore machines its uses snow to use two cores
  \item NOTE: Main bottleneck system memory: (per object: 3 * # of markers * # of phenotype * #of phenotype)
}
}
\value{
QCLscan object, a matrix with on the rows markers, and columns: phenotypes
}
\references{ 
  TODO
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  require(qtl)
  data(multitrait)
  multitrait = fill.geno(multitrait)
  res <- QCLscan(multitrait)
  plot(multitrait,res)
}

\keyword{methods}
