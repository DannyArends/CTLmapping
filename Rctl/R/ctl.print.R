#
# ctl.print.R
#
# copyright (c) 2010-2011 Danny Arends and Ritsert C. Jansen
# last modified Jan, 2012
# first written Oct, 2011
# 
# Print routines for CTL analysis
#

print.CTLobject <- function(x, ...){
  cat("CTLobject summary\n\n")
  cat("- Number of scanned phenotypes:",length(x),"\n")
}

getSignificantCTL <- function(CTLobject, threshold = 0.05, sep = ' '){
  all_significant <- vector("list", length(CTLobject))
  for(x in 1:length(CTLobject)){ #Get all significant CTLs
    significant_ctls <- names(which(apply(CTLobject[[x]]$ctl,1,function(x){any(x > -log10(threshold))})))
    if(length(significant_ctls) > 0){
      for(p in significant_ctls){
        cat(x, ctl.name(CTLobject[[x]]), p ,colnames(CTLobject[[x]]$ctl)[which.max(CTLobject[[x]]$ctl[p,])], max(CTLobject[[x]]$ctl[p,]),'\n',sep=sep)
      }
      all_significant[[x]] <- significant_ctls
    }
  }
  invisible(all_significant)
}

print.CTLscan <- function(x, ...){
  cat("CTLscan summary",attr(x$ctl,"name"),"\n\n")
  cat("- Number of background phenotypes",dim(x$ctl)[1],"\n")
  cat("- Number of markers",dim(x$ctl)[2],"\n")
  cat("- Number of permutations",length(unlist(x$p)),"\n")
}

getPermuteThresholds <-function(x, significance = c(.05,.01,.001), ..., verbose = FALSE){
  sorted <- sort(unlist(x))
  l <- length(sorted)
  values <- NULL
  valnames <- NULL
  for(x in significance){
    if(1/(x) < length(sorted)){
      v <- sorted[l*(1-(x))]
      values <- c(values,v)
      valnames <- c(valnames,paste(x*100,"%"))
      if(verbose) cat(x*100,"%\t",v,"\n")
    }else{
      values <- c(values,NaN)
      valnames <- c(valnames,paste(x*100,"%"))
      if(verbose) cat(x*100,"%\t",NaN,"\n")
    }
  }
  names(values) <- valnames
  values
}

print.CTLpermute <- function(x, ...){
  getPermuteThresholds(x, ..., verbose=TRUE)
}

# end of ctl.print.R
