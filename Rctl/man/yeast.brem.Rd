\name{yeast.brem}
\alias{yeast.brem}
\docType{data}
\title{ 
  Example gene expression data from Saccharomyces cerevisiae on 301 RNA expressions.
}
\description{
  Saccharomyces recombinant inbred lines. There are 109 lines, 301 phenotypes, genotyped 
  at 282 markers on 16 chromosomes stored as a list with 3 matrices: genotypes, phenotypes
  and map
}
\usage{
  data(yeast.brem)
}
\format{
  Data stored in a list holding 3 matrices genotypes, phenotypes and map
}
\details{
  Saccharomyces recombinant inbred lines. There are 109 lines, 301 RNA expression phenotypes.
  The individuals are genotyped at 282 markers on 16 chromosomes.
}
\source{
  Saccharomyces cerevisiae RNA expression data from XX, senior author: Rachel Brem 20XX, 
  Published in: Plos
}
\references{
  TODO
}
\examples{
  library(ctl)
  data(yeast.brem)              # Yeast data set
  
  yeast.brem$genotypes[1:5, 1:5]
  yeast.brem$phenotypes[1:5, 1:5]
  yeast.brem$map[1:5, ]
}
\keyword{datasets}
