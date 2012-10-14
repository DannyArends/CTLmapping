\name{CTLsignificant}
\alias{CTLsignificant}
\title{CTLsignificant - Get all significant interactions from a genome-wide CTLscan }
\description{
  Get all significant interactions from a genome-wide CTLscan.
}
\usage{
CTLsignificant <- function(CTLobject, significance = 0.05, what = c("names","ids"))
}

\arguments{
  \item{CTLobject}{ An object of class \code{"CTLobject"}, as output by \code{\link{CTLscan}}. }
  \item{significance}{ Significance threshold for setting a genomewide FDR. }
  \item{what}{ Return trait and marker names or column numbers (for indexing). }
}
\details{
  TODO
}
\value{
  A matrix significant CTL interactions with 4 columns: trait, marker, trait, lod
}
\references{ 
  TODO
}
\author{ 
  Danny Arends \email{Danny.Arends@gmail.com}\cr
  Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}
\note{
  TODO
}
\examples{
  library(ctl)
  data(ath.result)
  all_interactions <- CTLsignificant(ath.result)
  all_interactions[1:10, ]
  trait1_interactions <- CTLsignificant(ath.result[[1]])
  trait1_interactions
}

\keyword{methods}
