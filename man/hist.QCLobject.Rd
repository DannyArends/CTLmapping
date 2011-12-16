\name{hist.QCLobject}
\alias{hist.QCLobject}

\title{Plot histogram of QCL permutations}

\description{
  Plot histogram of QCL permutations (the output of \code{\link{QCLscan}}).
}

\usage{
\method{hist}{QCLobject}(x, pheno.col=1, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{pheno.col}{ Which phenotype column(s) should we analyse. Defaults to analyse all phenotype columns }
  \item{\dots}{ Passed to the function \code{\link[graphics]{image}} when it is called. }
}

\details{
  None.
}

\value{
  For a detailed description, see \code{\link{QCLprofiles}}
}

\examples{
  library(qcl)
  data(ath.metabolites) # Arabidopsis Thaliana dataset
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes, 1:5, n.perm=100)
  #Phenotype to phenotype matrix
  p2p_matrix <- image(qcl_result, against="phenotypes")
  #Phenotype to marker matrix
  p2m_matrix <- image(qcl_result, against="markers")
}

\seealso{
	\itemize{
    \item \code{\link{QCLscan}} - Scan for QCL
    \item \code{\link{QCLprofiles}} - Extract QCL interaction profiles
    \item \code{\link{print.QCLscan}} - Print a summary of a QCLscan
    \item \code{\link[graphics]{par}} - Plot parameters
    \item \code{\link[grDevices]{colors}} - Colors used in plotting
	}
}

\author{ 
	Danny Arends \email{Danny.Arends@gmail.com}\cr
	Maintainer: Danny Arends \email{Danny.Arends@gmail.com} 
}

\keyword{hplot}
