#
# ctl.statistics.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Jan, 2012
# first written Nov, 2011
# 
# R functions to do transform CTL mapping scores to Pvalues and LOD
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

CTLtoP <- function(CTLscan, onlySignificant = TRUE, verbose = TRUE){
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  permvalues <- sort(unlist(CTLscan$p))
  significant <- c()
  if(onlySignificant){
    maximums <- apply(abs(CTLscan$ctl),1,max)
    significant <- as.numeric(which(maximums > getPermuteThresholds(CTLscan$p)[1]))
  }  
  if(length(significant) > 1){
    scaled <- abs(CTLscan$ctl[significant, ])
    rnames <- rownames(CTLscan$ctl)[significant]
  }else{
    if(onlySignificant) warning("No significant found, converting all")
    scaled <- abs(CTLscan$ctl)
    rnames <- rownames(CTLscan$ctl)
  }
  pvalues <- unlist(lapply(1:length(permvalues),function(x){1-(x-1)/length(permvalues)}))
  result <- apply(scaled, 2, function(x){CTLtoPvalue.internal(x, permvalues, pvalues)})
  rownames(result) <- rnames
  result
}

CTLscoretoPvalue <- function(CTLscore, CTLpermute){
  if(missing(CTLscore)) stop("argument 'CTLscore' is missing, with no default")
  if(missing(CTLpermute)) stop("argument 'CTLpermute' is missing, with no default")
  permvalues <- sort(unlist(CTLpermute))
  pvalues <- unlist(lapply(1:length(permvalues),function(x){1-(x-1)/length(permvalues)}))
  CTLtoPvalue.internal(CTLscore, permvalues, pvalues)
}

#Determine a P-value based on the relative position of the score within the permutations
#Out of range values are tested using a GPD to estimate a P-value
CTLtoPvalue.internal <- function(CTLscore, permvalues, pvalues){
  cv <- as.numeric(permvalues[1])               #Critical value, every score < cv get set to p=1
  permvalues <- as.numeric(permvalues)
  ctlpvals <- unlist(lapply(CTLscore, function(score){
    score <- as.numeric(score)
    if(is.na(score)) return(1)                      #NA Value
    if(score < cv) return(1)                        #score < cv
    myrange <- which(permvalues > score)
    icx <- Inf
    if(length(myrange) > 0) icx <- min(myrange)
    if(is.finite(icx)){
      if(pvalues[min(icx)]==0) stop("P value of 0.0 at: pvalues[min(icx)]")
      return(pvalues[min(icx)])                     #Score falls inside the permutations
    }else{ # Score > any permutation p-value
      tryCatch(estimate <- extrapolateBeyondRange(permvalues, score),  
               error = function(e){ estimate <<- 1 })
      best_p <- pvalues[length(pvalues)]
      cat("Score ",score," is outside of permutation range: ",estimate,"/",best_p,"\n")
      if(estimate > best_p) return(best_p)
      return(estimate)
    }
  }))
  return(ctlpvals)
}

#Use the top 10% of permutation scores and fit a GPD uppertail distribution, then use the 
#GPD to obtain P-values for outliers, If the GPD estimates a 0 p-value we reduce our value 
#to get the minimum P-value
extrapolateBeyondRange <- function(permvalues, value = 0.6, top = 20){
  require(POT)
  gpd.threshold <- permvalues[(.80*length(permvalues))]
  mle <- fitgpd(permvalues, gpd.threshold, "mle")
  shape <- mle$param["shape"]
  scale <- mle$scale
  loc <- mle$threshold[1]
  dens <- function(x) dgpd(x, loc, scale, shape)
  warn <- FALSE
  prev.value <- value
  while(as.numeric(dens(value))==0 && value > 0){
    #cat("[FIX] Out of range and p=0:", value," ",as.numeric(dens(value)),"\n")
    warn <- TRUE
    value <- value - 0.0001
  }
  as.numeric(dens(value))
}

toLod <- function(CTLscan, onlySignificant = TRUE, verbose = TRUE){
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  ss <- proc.time()
  pm <- CTLtoP(CTLscan, onlySignificant, verbose)
  ee <- proc.time()
  if(verbose) cat("  - toLOD took",as.numeric(ee[3]-ss[3]),"seconds\n")
  -log10(pm)
}

CTLtoLODvector <- function(CTLscan, against = c("markers","phenotypes")){
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  if(!is.null(CTLscan$l)){
    if(against[1]=="markers")return(apply(CTLscan$l,2,sum))
    if(against[1]=="phenotypes")return(apply(CTLscan$l,1,max))
  }else{
    if(against[1]=="markers")return(apply(abs(CTLscan$ctl),2,sum))
    if(against[1]=="phenotypes")return(apply(abs(CTLscan$ctl),1,max))  
  }
}

CTLscantoScanone <- function(cross, CTLscan){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(missing(CTLscan)) stop("argument 'CTLscan' is missing, with no default")
  
  scores <- CTLtoLODvector(CTLscan)
  scores[which(!is.finite(scores))] <- NA
  lodscorestoscanone(cross, scores)
}

# end of ctl.statistics.R
