#
# ctl.power.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Pjotr Prins, Yang Li, and Ritsert C. Jansen
# last modified Dec, 2012
# first written Dec, 2012
# 
# Power analysis routines for CTL analysis
#

correlated.phenotype <- function(r = 0.5, nind = 5000, means = c(2,6), sds = c(0.3, 0.2)){
  sigma <- matrix(c(sds[1]^2, r*(sds[1]*sds[2]), r*(sds[1]*sds[2]), sds[2]^2), 2, 2)  
  invisible(mvrnorm(n=(nind*2), mu = means, Sigma = sigma))
}

check.marker <- function(marker, ratio = 50){
  if(sum(marker==1) <= 2) return(FALSE)
  if(sum(marker==2) <= 2) return(FALSE) 
  if( ((sum(marker==1) / sum(marker==2))/2) - (ratio/100) > 0.05) return(FALSE)
  return(TRUE)
}

check.pheno <- function(pheno, m0, m1, rdiff = 0.1, delta = 0.05){
  if(is.na(abs(cor(pheno[m0,1],pheno[m0,2]) - cor(pheno[m1,1],pheno[m1,2])-rdiff))) return(FALSE)
  if(abs(cor(pheno[m0,1],pheno[m0,2]) - cor(pheno[m1,1],pheno[m1,2])-rdiff) > delta) return(FALSE)
  return(TRUE)
}

create.marker <- function(nind, ratio = 50){
  marker <- as.numeric(runif(nind) > ratio / 100) + 1
  while(!check.marker(marker, ratio)){
    marker <- as.numeric(runif(nind) > ratio / 100) + 1
  }
  return(marker)
}

dcor.create <- function(rdiff = 0.1, nind = 5000, ratio = 50, delta = 0.05){
  marker <- create.marker(nind, ratio)
  m0 <- which(marker == 1);
  m1 <- which(marker == 2);
  pheno <- matrix(0, nind, 2)
  while(!check.pheno(pheno, m0, m1, rdiff, delta)){
    r0 <- runif(1)
    r1 <- r0-rdiff
    while(r1 < -1){
      r0 <- runif(1)
      r1 <- r0-rdiff
    }
    p0 <- correlated.phenotype(r0, length(m0)/2)
    p1 <- correlated.phenotype(r1, length(m1)/2)
    pheno[m0,1] <- p0[,1]
    pheno[m0,2] <- p0[,2]

    pheno[m1,1] <- p1[,1]
    pheno[m1,2] <- p1[,2]
  }
  return(list(pheno, t(t(marker))))
}

getRatio <- function(majorFreq = 0.9){
  minFreq = 1-majorFreq
  pp <- majorFreq * majorFreq
  pq <- 2 * minFreq * majorFreq
  qq <- minFreq * minFreq
  100*((pq+qq)/pp)
}

CTLpowerAnalysis <- function(n, effects=seq(0, 1, 0.05), individuals=c(20, 40, 60), ratios=seq(10, 50, 5), ...){
  output <- NULL
  for(rat in ratios){  for(eff in effects){  for(ind in individuals){
    sign_cnt = 0
    for(x in 1:n){
      input   <- dcor.create(eff, ind, rat)
      ctl_res <- CTLscan(input[[2]], input[[1]], ... , verbose=FALSE)
      if(!is.nan(ctl_res[[1]]$ctl[2])){
        if(ctl_res[[1]]$ctl[2] > -log10(0.05)){ 
          sign_cnt = sign_cnt + 1; 
        }
      }
      #cat("Done ",x,"\n")
    }
    output <- rbind(output, c(rat, eff, ind, sign_cnt / n))
    cat("Done with:",rat, "/", eff, "on", ind, "Individuals\n")
  }}}
  return(output)
}

CTLanalyseEffects <- function(n=1000, perms=100, effects=seq(0, 1, 0.05), ind = 100, ratio=50, ...){
  output <- NULL
  for(eff in effects){
    catout <- NULL
    for(p in 1:perms){
      sign_cnt = 0
      for(x in 1:n){
        input   <- dcor.create(eff, ind, ratio)
        ctl_res <- CTLscan(input[[2]], input[[1]], ... , verbose=FALSE)
        if(!is.nan(ctl_res[[1]]$ctl[2])){
          if(ctl_res[[1]]$ctl[2] > -log10(0.05)){ 
            sign_cnt = sign_cnt + 1; 
          }
        }
      }
      catout <- c(catout, (sign_cnt/n))
      cat("Done ",p,"/",perms,"\n")
    }
    output <- cbind(output, catout)
  }
  colnames(output) <- effects
  return(output)
}

CTLanalyseRatios <- function(n = 100, perms = 100, effect = 0.4, ind = 200, ratios=c(5, 10, 20, 30, 40, 50), ...){
  output <- NULL
  for(ratio in ratios){
    catout <- NULL
    for(p in 1:perms){
      sign_cnt = 0
      for(x in 1:n){
        input   <- dcor.create(effect, ind, ratio)
        ctl_res <- CTLscan(input[[2]], input[[1]], ... , verbose=FALSE)
        if(!is.nan(ctl_res[[1]]$ctl[2])){
          if(ctl_res[[1]]$ctl[2] > -log10(0.05)){ 
            sign_cnt = sign_cnt + 1; 
          }
        }
      }
      catout <- c(catout, (sign_cnt/n))
      cat("Done ",p,"/",perms,"\n")
    }
    output <- cbind(output, catout)
  }
  colnames(output) <- ratios
  return(output)
}

test.power.test <- function(){
  sample100 <- CTLanalyseEffects(effects=c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75), ind=100)
  sample200 <- CTLanalyseEffects(effects=c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75), ind=200)
  sample400 <- CTLanalyseEffects(effects=c(0.02, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75), ind=400)
  
  sample100r <- CTLanalyseRatios(ind=100)
  sample200r <- CTLanalyseRatios(ind=200)
  sample400r <- CTLanalyseRatios(ind=400)

  setwd("~/Github/Articles/CTLpaper/nar/img")
  postscript("power_analysis.eps",width=12, height=6,paper="special", horizontal=FALSE)
  op <- par(mfrow=c(1,2))

  plot(c(0, 0.8),c(0, 1),t='n',xlab="Simulated effect size",ylab="Power", xaxt="n")
  boxplot(sample100, at=as.numeric(colnames(sample100)), boxwex=0.02, add=T,pch=20)
  points(as.numeric(colnames(sample100)), apply(sample100,2,median), t='l')
  boxplot(sample200, at=as.numeric(colnames(sample200)), boxwex=0.02, add=T, col="blue",pch=20)
  points(as.numeric(colnames(sample200)), apply(sample200,2,median), t='l', col="blue")
  boxplot(sample400, at=as.numeric(colnames(sample400)), boxwex=0.02, add=T, col="orange",pch=20)
  points(as.numeric(colnames(sample400)), apply(sample400,2,median), col='orange',t='l')

  legend("bottomright",c("100","200","300"),col=c("black","blue","orange"),lwd=1,title="Sample size")

  plot(c(0, 50),c(0, 1),t='n',xlab="Minor allele frequency",ylab="Power", xaxt="n")
  boxplot(sample100r, at=as.numeric(colnames(sample100r)), boxwex=1, add=T,pch=20)
  points(as.numeric(colnames(sample100r)), apply(sample100r,2,median), t='l')
  boxplot(sample200r, at=as.numeric(colnames(sample200r)), boxwex=1, add=T, col="blue",pch=20)
  points(as.numeric(colnames(sample200r)), apply(sample200r,2,median), t='l', col="blue")
  boxplot(sample400r, at=as.numeric(colnames(sample400r)), boxwex=1, add=T, col="orange",pch=20)
  points(as.numeric(colnames(sample400r)), apply(sample400r,2,median), col='orange',t='l')

  legend("bottomright",c("100","200","300"),col=c("black","blue","orange"),lwd=1,title="Sample size")
  dev.off()
}

# end of ctl.power.R
