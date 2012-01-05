#
# qcl.hist.R
#
# copyright (c) 2010-2011 Danny Arends and Ritsert C. Jansen
# last modified Dec, 2011
# first written Dec, 2011
# 
# Histogram routines for QCL analysis
#

hist.QCLobject <- function(x, pheno.col=1, ...){
  if(missing(x)) stop("argument 'x' which expects a 'QCLobject' object is missing, with no default")
  namez <- NULL
  plot(c(0,0.8),c(0,length(unlist(x[[pheno.col[1]]]$p))/5),type='n',main="QCL scores permutations",ylab="Frequency",xlab="Difference in correlation^2")
  for(pheno in pheno.col){
    if(is.null(x[[pheno]]$p)) stop(paste("Permutations not found for pheno.col=",pheno))
    sorted <- sort(unlist(x[[pheno]]$p))
    hist(sorted, breaks=seq(0,0.8,0.01), add=TRUE, col=pheno, ...)
    namez <- c(namez,attr(x[[pheno]]$qcl,"name"))
  }
  legend("topright", namez, col=pheno.col,lwd=6)
}
