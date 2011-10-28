\name{download.KEGG.annotation}
\alias{download.KEGG.annotation}
\title{download.KEGG.annotation - Download gene and pathway information from KEGG }
\description{
	Download gene and pathway information from KEGG
}
\usage{
download.KEGG.annotation(genenames, species = "sce", directory = "tmp_kegg", force.download=FALSE, verbose = FALSE)
}

\arguments{
  \item{genenames}{ A vector holding the names of the genes to be annotated }
  \item{species}{ KEGG Species code, see: \href{http://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=mode.map}{KEGG organisms} }
  \item{directory}{ Directory to store the downloaded files (used as cache) }
  \item{force.download}{ Force the download of all files, overwrites the cached files }
  \item{verbose}{ Be verbose }
}
\details{
Default behaviour is force.download = FALSE, which creates a local cache used for annotation
\itemize{
  \item NOTE: NO check for unknown names, just searches
}
}
\value{
  List of length genenames with at each index the KEGG pathway annotation of the gene
}
\references{
  \href{http://www.genome.jp/kegg/kegg2.html}{KEGG database}
}
\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}\cr
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  pathwaydata <- download.KEGG.annotation(c("YLR075W","YER110C","YOR210W"),verbose=TRUE)
}

\keyword{methods}
