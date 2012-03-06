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
  l <- length(permvalues)
  if(onlySignificant){
    mysignificant <- as.numeric(which(apply(abs(CTLscan$ctl),1,max) > getPermuteThresholds(CTLscan$p)[1]))
    if(length(mysignificant) > 1){
      scaled <- abs(CTLscan$ctl[mysignificant, ])
      rnames <- rownames(CTLscan$ctl)[mysignificant]
    }else{
      scaled <- abs(CTLscan$ctl)
      rnames <- rownames(CTLscan$ctl)
    }
  }else{
    scaled <- abs(CTLscan$ctl)
    rnames <- rownames(CTLscan$ctl)
  }
  result <- apply(scaled, 2, function(x){CTLtoPvalue.internal(x, permvalues, l)})
  rownames(result) <- rnames
  result
}

CTLscoretoPvalue <- function(CTLscore, CTLpermute){
  if(missing(CTLscore)) stop("argument 'CTLscore' is missing, with no default")
  if(missing(CTLpermute)) stop("argument 'CTLpermute' is missing, with no default")
  permvalues <- sort(unlist(CTLpermute))
  l <- length(permvalues)
  CTLtoPvalue.internal(CTLscore, permvalues, l)
}

#Determine a P-value based on the relative position of the score within the permutations
#Out of range values are tested using a GPD to estimate a P-value
CTLtoPvalue.internal <- function(CTLscore, permvalues, l){
  res <- NULL
  warn <- TRUE
  for(y in CTLscore){
    if(!is.na(which(permvalues > y)&&1)){
      index_l <- min(which(permvalues > y))
    }else{
      index_l <- Inf
    }
    if(!is.na(which(permvalues < y)&&1)){
      index_r <- max(which(permvalues < y))
    }else{
      index_r <- Inf
    }
    if(!is.finite(index_l)){
      tryCatch(estimate <- extrapolateBeyondRange(permvalues, y),  error = function(e) {estimate <<- 1})
      if(estimate > 1-((l-1)/l)){
        estimate <- 1-((l-1)/l)
      }
      res <- c(res,estimate)
      if(warn){
        cat("  - [Warning] Scores out of permutation range, please do more permutations\n")
        warn <- FALSE  
      }
    }
    if(!is.finite(index_r)){
     estimate <- (1-(index_l/l))
      if(estimate==0){
        res <- c(res,1)
      }else{
        res <- c(res,estimate)
      }
    }
    if(is.finite(index_l) && is.finite(index_r)){
      res <- c(res,1-(index_r/l))
    }
  }
  res
}

#Use the top 10% of permutation scores and fit a GPD uppertail distribution, then use the 
#GPD to obtain P-values for outliers, If the GPD is 0 we reduce our value till we get the 
#minimum P-value
extrapolateBeyondRange <- function(permvalues, value = 0.6){
  require(POT)
  gpd.threshold <- permvalues[(.80*length(permvalues))]
  mle <- fitgpd(permvalues, gpd.threshold, "mle")
  shape <- mle$param["shape"]
  scale <- mle$scale
  loc <- mle$threshold[1]
  dens <- function(x) dgpd(x, loc, scale, shape)
  warn <- FALSE
  prev.value <- value
  while(as.numeric(dens(value))==0){
    warn <- TRUE
    value <- value - 0.0001
  }
  #if(warn){
  ##  cat("Warning: scores out of permutation range, unable to estimate correctly",value,"/",prev.value,dens(value),"\n")
  #}
  as.numeric(dens(value))
}

toLod <- function(CTLscan, onlySignificant = TRUE, verbose = FALSE){
  ss <- proc.time()
  pmatrix <- CTLtoP(CTLscan, onlySignificant, verbose)
  ee <- proc.time()
  if(verbose) cat("  - toLOD took",as.numeric(ee[3]-ss[3]),"seconds\n")
  -log10(pmatrix)
}

CTLtoLODvector <- function(CTLscan, against = c("markers","phenotypes")){
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
