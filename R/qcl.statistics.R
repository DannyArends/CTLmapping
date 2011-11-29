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

QCLtoPvalue <- function(QCLscan, permutations, pheno.col=1, onlySignificant = TRUE){
  if(missing(permutations)) stop("Please provide permutations")
  maximums <- lapply(permutations[[pheno.col]], function(x){
    apply(abs(x), 2, max)
  })
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  if(onlySignificant){
    mysignificant <- as.numeric(which(apply(abs(QCLscan[[pheno.col]]),1,max) > getPermuteThresholds(QCL_p, pheno.col)[1]))
    if(length(mysignificant) != 0){
      scaled <- abs(QCLscan[[pheno.col]][mysignificant, ])
    }else{
      scaled <- abs(QCLscan[[pheno.col]])
    }
  }else{
    scaled <- abs(QCLscan[[pheno.col]])
  }
  apply(scaled, 1, function(x){QCLtoPvalue.internal(x,sorted,l)})
}

QCLscoretoPvalue <- function(QCLscore, permutations){
  if(missing(permutations)) stop("Please provide permutations")
  maximums <- lapply(permutations, function(x){
    apply(abs(x), 2, max)
  })
  sorted <- sort(unlist(maximums))
  l <- length(sorted)
  QCLtoPvalue.internal(QCLscore,sorted,l)
}

QCLtoPvalue.internal <- function(QCLscore, sorted, l){
  res <- NULL
  warn <- TRUE
  for(y in QCLscore){
    if(!is.na(which(sorted > y)&&1)){
      index_l <- min(which(sorted > y))
    }else{
      index_l <- Inf
    }
    if(!is.na(which(sorted < y)&&1)){
      index_r <- max(which(sorted < y))
    }else{
      index_r <- Inf
    }
    if(!is.finite(index_l)){
      res <- c(res,extrapolateBeyondRange(sorted, y)) #res <- c(res,1-((index_r-1)/l))
      if(warn){
        cat("Warning: scores out of permutation range, please do more permutations\n")
        warn <- FALSE  
      }
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

extrapolateBeyondRange <- function(permvalues, value = 0.6){
  require(POT)
  #We use the top 10% of data to estimate the GPD uppertail distribution
  mle <- fitgpd(permvalues, permvalues[.90*length(permvalues)], "mle")
  shape <- mle$param["shape"]
  scale <- mle$scale
  loc <- mle$threshold[1]
  dens <- function(x) dgpd(x, loc, scale, shape)
  warn <- FALSE
  prev.value <- value
  while(as.numeric(dens(value))==0){
    warn <- TRUE
    value <- value - 0.01
  }
  if(warn){
    cat("Warning: scores out of permutation range, unable to estimate correctly",value,"/",prev.value,"\n")
  }
  dens(value)
}

QCLtoLOD <- function(QCLscan, permutations, pheno.col = 1, onlySignificant = TRUE){
  -log10(QCLtoPvalue(QCLscan, permutations, pheno.col, onlySignificant))
}

QCLtoLODvector <- function(QCLscan, permutations, pheno.col = 1){
  apply(QCLtoLOD(QCLscan, permutations, pheno.col, FALSE),1,sum)
}

QCLscantoScanone <- function(cross, QCLscan, permutations, pheno.col=1){
  scores <- QCLtoLODvector(QCLscan, permutations, pheno.col)
  scores[which(!is.finite(scores))] <- NA
  lodscorestoscanone(cross, scores)
}

plot.QCLpermute <- function(x, ...){
  plot(seq(0,0.5,0.01),QCLscoretoPvalue(seq(0,0.5,0.01),x),main="QCL to P.value",xlab="QCL",ylab="Pvalue")
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
  if(missing(x)) stop("Please provide permutations")
  maximums <- lapply(x, function(v){
    apply(abs(v), 2, max)
  })
  sorted <- sort(unlist(maximums))
  hist(sorted,breaks=100,main="QCL scores during permutation",...)
}
