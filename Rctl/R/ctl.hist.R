#
# ctl.hist.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Dec, 2011
# 
# Histogram routines for CTL analysis

hist.CTLobject <- function(x, pheno.col=1, ...){
  if(missing(x)) stop("argument 'x' which expects a 'CTLobject' object is missing, with no default")
  namez <- NULL
  plot(c(0,1.0),c(0,length(unlist(x[[pheno.col[1]]]$perms))/5),type='n',main="CTL scores permutations",ylab="Frequency",xlab="Difference in correlation^2")
  for(pheno in pheno.col){
    if(is.null(x[[pheno]]$perms)) stop(paste("Permutations not found for pheno.col=",pheno))
    sorted <- sort(unlist(x[[pheno]]$perms))
    hist(sorted, breaks=seq(0,1.0,0.01), add=TRUE, col=pheno, ...)
    namez <- c(namez,ctl.name(x[[pheno]]))
  }
  legend("topright", legend=namez, col=pheno.col,lwd=6)
}

# end of ctl.hist.R
