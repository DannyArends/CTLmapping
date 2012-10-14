\name{ath.churchill}
\alias{ath.churchill}
\alias{ath.gary}
\docType{data}

\title{ 
  Example metabolite expression data from Arabidopsis Thaliana on 9 metabolites.
}

\description{
  Arabidopsis recombinant inbred lines by selfing.  There are 403 lines, 9 phenotypes, 
  and 69 markers on 5 chromosomes stored as a list with 3 matrices: genotypes, phenotypes, map
}
\usage{
  data(ath.churchill)
}
\format{
  Data stored in a list holding 3 matrices genotypes, phenotypes and map
}
\details{
  Arabidopsis recombinant inbred lines by selfing.  There are 403 lines, 9 metabolic 
  phenotypes, and 69 markers on 5 chromosomes.
}
\source{
  Arabidopsis Bay-0 x Sha metabolite data from XX, senior author: Gary Churchill 2012, 
  Published in: Plos
}
\references{
  TODO
}
\examples{
  library(ctl)
  data(ath.churchill)           # Arabidopsis thaliana dataset
  
  ath.gary$genotypes[1:5, 1:5]  # ath.gary is the short name
  ath.gary$phenotypes[1:5, 1:5]
  ath.gary$map[1:5, ]
}
\keyword{datasets}
