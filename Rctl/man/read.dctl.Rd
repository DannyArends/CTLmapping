\name{read.dctl}
\alias{read.dctl}

\title{Print the results of a CTL genome scans}

\description{
  Print the results of a \code{\link{CTLscan}}
}

\usage{
read.dctl(geno.file="genotypes.csv", pheno.file="phenotypes.csv", results="output/", verbose = TRUE)
}

\arguments{
  \item{geno.file}{ Matrix of genotypes. (markers x individuals) }
  \item{pheno.file}{ Matrix of phenotypes. (phenotypes x individuals) }
  \item{results}{ Directory in which CTL results from D are located. }
  \item{verbose}{ Be verbose }
}

\details{
  None.
}

\value{
  None. 
}

\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\examples{
  setwd("E:/github/Rpackages/CTLmapping/D/")
  ctl_res <- read.dctl("test/data/genotypes.csv", "test/data/phenotypes.csv", "output/")
  image(ctl_res)
}

\seealso{
  \itemize{
    \item \code{\link{CTLscan}} - Scan for CTL
    \item \code{\link{plot.CTLscan}} - Plot a CTLscan object
    \item \code{\link{image.CTLobject}} - Plot a CTLscan object
  }
}

\references{
  TODO
}

\keyword{print}
