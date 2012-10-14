#
# ctl.print.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Oct, 2012
# first written Oct, 2011
# 
# Print routines for CTL analysis
#

print.CTLobject <- function(x, ...){
  cat("CTLobject summary\n\n")
  cat("- Number of scanned phenotypes:",length(x),"\n")
}

CTLsignificant <- function(CTLobject, significance = 0.05, what = c("names","ids")){
  if(any(class(CTLobject)=="CTLscan")) CTLobject = list(CTLobject)
  all_sign <- NULL
  if(length(what) > 1) what = what[1]
  for(x in 1:length(CTLobject)){ #Get all significant CTLs
    p_above <- which(apply(CTLobject[[x]]$ctl,2,function(x){any(x > -log10(significance))}))
    if(what != "ids"){ p_above <- names(p_above) }
    
    if(length(p_above) > 0){
      for(p in p_above){
        m_above <- which(CTLobject[[x]]$ctl[,p] > -log10(significance))
        if(what != "ids"){ m_above <- names(m_above) }
        for(m in m_above){
          if(what == "ids"){
            all_sign <- rbind(all_sign, c(x, m, p, CTLobject[[x]]$ctl[m, p]) )
          }else{
            all_sign <- rbind(all_sign, c(ctl.name(CTLobject[[x]]), m, p, CTLobject[[x]]$ctl[m, p]) )
          }
        }
      }
    }
  }
  items <- 0
  if(!is.null(all_sign)){
    all_sign <- as.data.frame(all_sign)
    all_sign[,4] <- round(as.numeric(as.character(all_sign[,4])),d=2)
    colnames(all_sign) <- c("trait","marker","trait","lod")
    items <- nrow(all_sign)
  }
  cat("Found",items,"significant CTLs\n")
  invisible(all_sign)
}

print.CTLscan <- function(x, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
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
