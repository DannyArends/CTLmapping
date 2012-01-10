#
# qcl.statistics.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Jan, 2012
# first written Nov, 2011
# 
# R functions to do transform QCL mapping scores to Pvalues and LOD
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

QCLtoP <- function(QCLscan, onlySignificant = TRUE, verbose = TRUE){
  if(missing(QCLscan)) stop("argument 'QCLscan' is missing, with no default")
  permvalues <- sort(unlist(QCLscan$p))
  l <- length(permvalues)
  if(onlySignificant){
    mysignificant <- as.numeric(which(apply(abs(QCLscan$qcl),1,max) > getPermuteThresholds(QCLscan$p)[1]))
    if(length(mysignificant) > 1){
      scaled <- abs(QCLscan$qcl[mysignificant, ])
      rnames <- rownames(QCLscan$qcl)[mysignificant]
    }else{
      scaled <- abs(QCLscan$qcl)
      rnames <- rownames(QCLscan$qcl)
    }
  }else{
    scaled <- abs(QCLscan$qcl)
    rnames <- rownames(QCLscan$qcl)
  }
  result <- apply(scaled, 2, function(x){QCLtoPvalue.internal(x, permvalues, l)})
  rownames(result) <- rnames
  result
}

QCLscoretoPvalue <- function(QCLscore, QCLpermute){
  if(missing(QCLscore)) stop("argument 'QCLscore' is missing, with no default")
  if(missing(QCLpermute)) stop("argument 'QCLpermute' is missing, with no default")
  permvalues <- sort(unlist(QCLpermute))
  l <- length(permvalues)
  QCLtoPvalue.internal(QCLscore, permvalues, l)
}

#Determine a P-value based on the relative position of the score within the permutations
#Out of range values are tested using a GPD to estimate a P-value
QCLtoPvalue.internal <- function(QCLscore, permvalues, l){
  res <- NULL
  warn <- TRUE
  for(y in QCLscore){
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

toLod <- function(QCLscan, onlySignificant = TRUE, verbose = FALSE){
  ss <- proc.time()
  pmatrix <- QCLtoP(QCLscan, onlySignificant, verbose)
  ee <- proc.time()
  if(verbose) cat("  - toLOD took",as.numeric(ee[3]-ss[3]),"seconds\n")
  -log10(pmatrix)
}

QCLtoLODvector <- function(QCLscan, against = c("markers","phenotypes")){
  if(!is.null(QCLscan$l)){
    if(against[1]=="markers")return(apply(QCLscan$l,2,sum))
    if(against[1]=="phenotypes")return(apply(QCLscan$l,1,sum))
  }else{
    if(against[1]=="markers")return(apply(abs(QCLscan$qcl),2,sum))
    if(against[1]=="phenotypes")return(apply(abs(QCLscan$qcl),1,sum))  
  }
}

QCLscantoScanone <- function(cross, QCLscan){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(missing(QCLscan)) stop("argument 'QCLscan' is missing, with no default")
  
  scores <- QCLtoLODvector(QCLscan)
  scores[which(!is.finite(scores))] <- NA
  lodscorestoscanone(cross, scores)
}

# end of qcl.statistics.R
