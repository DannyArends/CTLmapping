  
correlated.phenotype <- function(r = 0.5, nind = 5000, means = c(2,6), sds = c(0.3, 0.2)){
  a <- rmvnorm(n=(nind*2), mu=means, sigma = matrix(c(sds[1]^2,r*(sds[1]*sds[2]),r*(sds[1]*sds[2]),sds[2]^2), 2, 2)) 
  invisible(a)
}

create.marker <- function(nind, ratio = 50){
  marker <- as.numeric(runif(nind) > ratio / 100)+1
  while(sum(marker==1) <= 2 || sum(marker==2) <= 2 || ((sum(marker==1)/sum(marker==2))/2) - (ratio/100) > 0.05){
    marker <- as.numeric(runif(nind) > ratio / 100)+1
  }
  return(marker)
}

dcor.create <- function(rdiff = 0.1, nind = 5000, ratio = 50, delta = 0.05){
  marker <- create.marker(nind, ratio)
  m0 <- which(marker == 1);
  m1 <- which(marker == 2);
  pheno <- matrix(0, nind, 2)
  while(is.na(abs(cor(pheno[m0,1],pheno[m0,2]) - cor(pheno[m1,1],pheno[m1,2])-rdiff)) || abs(cor(pheno[m0,1],pheno[m0,2]) - cor(pheno[m1,1],pheno[m1,2])-rdiff) > delta){
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
  #cat(abs(cor(pheno[m0,1],pheno[m0,2]) - cor(pheno[m1,1],pheno[m1,2])), "\n")
  return(list(pheno, marker))
}

power.test <- function(n = 1000, effects = seq(0, 1, 0.05), individuals = c(20, 40, 60), ratios = seq(10, 50, 5), method = c("Exact","Power","Adjacency")){
  output <- NULL  
  for(rat in ratios){  
  for(eff in effects){
  for(ind in individuals){
    sign_cnt = 0
    for(x in 1:n){
      input   <- dcor.create(eff, ind, rat) 
      ctl_res <- CTLscan(t(t(input[[2]])), input[[1]], method=method[1], n.perms = 500, n.cores=1, verbose=FALSE)
      if(ctl_res[[1]]$ctl[2] > -log10(0.05)){ sign_cnt = sign_cnt + 1; }
    }
    output <- rbind(output, c(rat, eff, ind, sign_cnt / n))
    cat("Done with:",rat, "/", eff, "on",ind,"Individuals\n")
  } } }
  return(output)
}

test.power.test <- function(){
  library(ctl)
  library(mvtnorm)
  res1 <- power.test(100, individuals = c(100),method="Exact")
  res2 <- power.test(100, individuals = c(100),method="Power")
  res3 <- power.test(100, individuals = c(100),method="Adjacency")

  plot(c(5,55),c(0,1),t='n', ylab="Effect Size", xlab="Genotype ratio")
  xdist <- 5; ydist <- 0.05
  for(x in 1:nrow(res1)){ 
    rect(res1[x,1]-(xdist/2), res1[x,2]-(ydist/2), res1[x,1]+(xdist/2), res1[x,2]+(ydist/2),col=rgb(1-res1[x,4], res1[x,4], 0), border="white"); 
  }
  box()
  plot(c(5,55),c(0,1),t='n', ylab="Effect Size", xlab="Genotype ratio")
  xdist <- 5; ydist <- 0.05
  for(x in 1:nrow(res2)){ 
    rect(res2[x,1]-(xdist/2), res2[x,2]-(ydist/2), res2[x,1]+(xdist/2), res2[x,2]+(ydist/2),col=rgb(1-res2[x,4], res2[x,4], 0), border="white"); 
  }
  box()
  plot(c(5,55),c(0,1),t='n', ylab="Effect Size", xlab="Genotype ratio")
  xdist <- 5; ydist <- 0.05
  for(x in 1:nrow(res3)){ 
    rect(res3[x,1]-(xdist/2), res3[x,2]-(ydist/2), res3[x,1]+(xdist/2), res3[x,2]+(ydist/2),col=rgb(1-res3[x,4], res3[x,4], 0), border="white"); 
  }
  box()

}
