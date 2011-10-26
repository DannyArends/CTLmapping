\name{download.KEGG.annotation}
\alias{download.KEGG.annotation}
\title{download.KEGG.annotation - Download pathway information from KEGG }
\description{
	Download pathway information from KEGG
}
\usage{
download.KEGG.annotation(genenames, species = "sce", directory = "tmp_kegg", force.download=FALSE, verbose = FALSE)
}

\arguments{
  \item{genenames}{ A vector holding the names of the genes to be annotated }
  \item{species}{ Species specific gene code, see: \url{http://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=mode.map}{KEGG organisms} }
  \item{directory}{ Directory to store the downloaded files (used as cache) }
  \item{force.download}{ Force downloading of files, do not use the cached versions }
  \item{verbose}{ Be verbose }
}
\details{
TODO
\itemize{
  \item NOTE: NO check for unknown names, just searches
}
}
\value{
  List of length genenames with at each index the KEGG pathway annotation of the gene
}
\references{
  \url{http://www.genome.jp/kegg/kegg2.html}{KEGG database}
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  pathwaydata <- download.KEGG.annotation(c("YLR075W","YER110C","YOR210W"),verbose=TRUE)
}

\keyword{methods}
