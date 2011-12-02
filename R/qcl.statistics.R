#
# qcl.statistics.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Nov, 2011
# first written Nov, 2011
# 
# R functions to do transform QCL mapping scores to Pvalues and LOD
# Example data C. Elegans and available at request ( Danny.Arends@gmail.com )
#

QCLtoPvalue <- function(QCLscan, QCLpermute, pheno.col=1, onlySignificant = TRUE){
  if(missing(QCLscan)) stop("argument 'QCLscan' is missing, with no default")
  if(missing(QCLpermute)) stop("argument 'QCLpermute' is missing, with no default")
  permvalues <- sort(unlist(QCLpermute[[pheno.col]]))
  l <- length(permvalues)
  if(onlySignificant){
    mysignificant <- as.numeric(which(apply(abs(QCLscan[[pheno.col]]),1,max) > getPermuteThresholds(QCLpermute, pheno.col)[1]))
    if(length(mysignificant) > 1){
      scaled <- abs(QCLscan[[pheno.col]][mysignificant, ])
      rnames <- rownames(QCLscan[[pheno.col]])[mysignificant]
    }else{
      scaled <- abs(QCLscan[[pheno.col]])
      rnames <- rownames(QCLscan[[pheno.col]])
    }
  }else{
    scaled <- abs(QCLscan[[pheno.col]])
    rnames <- rownames(QCLscan[[pheno.col]])
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
      estimate <- extrapolateBeyondRange(permvalues, y)
      if(estimate > 1-((l-1)/l)){
        estimate <- 1-((l-1)/l)
      }
     # cat(y,estimate,"\n")
      res <- c(res,estimate)
     # if(warn){
     #   cat("Warning: scores out of permutation range, please do more permutations\n")
     #   warn <- FALSE  
     # }
    }
    if(!is.finite(index_r)){
      res <- c(res,1-(index_l/l))
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

QCLtoLOD <- function(QCLscan, QCLpermute, pheno.col = 1, onlySignificant = TRUE){
  -log10(QCLtoPvalue(QCLscan, QCLpermute, pheno.col, onlySignificant))
}

QCLtoLODvector <- function(QCLscan, QCLpermute, pheno.col = 1, against = c("markers","phenotypes")){
  if(against[1]=="markers")return(apply(QCLtoLOD(QCLscan, QCLpermute, pheno.col, FALSE),2,sum))
  if(against[1]=="phenotypes")return(apply(QCLtoLOD(QCLscan, QCLpermute, pheno.col, FALSE),1,sum))
}

QCLscantoScanone <- function(cross, QCLscan, QCLpermute, pheno.col=1){
  if(missing(cross)) stop("argument 'cross' is missing, with no default")
  if(missing(QCLscan)) stop("argument 'QCLscan' is missing, with no default")
  if(missing(QCLpermute)) stop("argument 'QCLpermute' is missing, with no default")
  
  scores <- QCLtoLODvector(QCLscan, QCLpermute, pheno.col)
  scores[which(!is.finite(scores))] <- NA
  lodscorestoscanone(cross, scores)
}

plot.QCLpermute <- function(x, ...){
  if(missing(x)) stop("argument 'x' which expects a 'QCLpermute' object is missing, with no default")
  plot(seq(0,0.9,0.01),QCLscoretoPvalue(seq(0,0.9,0.01),x),main="QCL to P.value",xlab="QCL",ylab="Pvalue")
  significant <- print.QCLpermute(x)
  mycolors <- c("red","orange","green")
  idx <- 1
  for(y in c(.05,.01,.001)){
    lines(rbind(c(-1,y),c(significant[idx],y)),lty=2,col=mycolors[idx])
    lines(rbind(c(significant[idx],-1),c(significant[idx],y)),lty=2,col=mycolors[idx])
    idx <- idx+1
  }
}

hist.QCLpermute <- function(x, ...){
  if(missing(x)) stop("argument 'x' which expects a 'QCLpermute' object is missing, with no default")
  maximums <- lapply(x, function(v){
    apply(abs(v), 2, max)
  })
  sorted <- sort(unlist(maximums))
  hist(sorted,breaks=100,main="QCL scores during permutation",...)
}
