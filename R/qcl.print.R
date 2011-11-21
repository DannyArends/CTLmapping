#
# qcl.print.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Oct, 2011
# 
# Plotting routines for QCL analysis
#

print.QCLscan <- function(x, ...){
  cat("QCLscan summary\n\n")
  cat("- Number of phenotypes scanned",length(x),"/",dim(x[[1]])[1],"\n")
  cat("- Number of markers",dim(x[[1]])[2],"\n")
}

print.QCL  <- function(x, ...){
  cat("Phenotype",attr(x,"name"),"QCL summary\n")
  highest <- 0
  np <- NULL
  nm <- NULL
  cutoffs <- c(0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65)
  for(t in cutoffs){
    smprofile <- QCLprofile(x, t, FALSE)
    if(!is.null(dim(smprofile)[1])) highest = t
    np  <- c(np,dim(smprofile)[1])
    nm  <- c(nm,dim(smprofile)[2])
  }
  cat("Threshold:",  cutoffs,"\n",sep="\t")
  cat("Phenotypes:", np,"\n",sep="\t")
  cat("Markers:",    nm,"\n\n",sep="\t")
  cat("Interaction matrix at highest cut-off", highest,":\n")
  QCLprofile(x, highest, TRUE)
}

print.QCLpermute <- function(x, ...){
  maximums <- lapply(x,function(v){
    apply(abs(v),2,max)}
  )
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in c(.95,.99,.999)){
    values <- c(values,sorted[l*x])
    valnames <- c(valnames,paste((1-x)*100,"%"))
    cat((1-x)*100,"%\t",sorted[l*x],"\n")
  }
  names(values) <- valnames
  invisible(values)
}
