\name{image.QCLscan}
\alias{image.QCLscan}

\title{Plot genome-wide QCL on multiple traits}

\description{
  Plot the QCLs for genome-wide QCL on multiple traits (the output of \code{\link{QCLscan}}).
}

\usage{
\method{image}{QCLscan}(x, QCLpermute, against = c("markers","phenotypes"), 
        do.grid = TRUE, grid.col = "black", verbose=FALSE, \dots)
}

\arguments{
  \item{x}{ An object of class \code{"QCLscan"}, as output by \code{\link{QCLscan}}. }
  \item{QCLpermute}{ An object of class \code{"QCLpermute"}, as output by \code{\link{QCLpermute}}. }
  \item{against}{ Image one of these interaction matrices, 
                  this parameter is passed to the \code{\link{QCLprofiles}} function }
  \item{do.grid}{ Draw grid lines. }
  \item{grid.col}{ Colors used to draw the grid. }
    \item{verbose}{ Be verbose }
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
  qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes)
  qcl_perms <- QCLpermute(ath.metab$genotypes, ath.metab$phenotypes,1:24,15)
  #Phenotype to phenotype matrix
  p2p_matrix <- image(qcl_result, qcl_perms, against="phenotypes")
  #Phenotype to marker matrix
  p2m_matrix <- image(qcl_result, qcl_perms, against="markers")
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
