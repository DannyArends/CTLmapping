tolist <- function(x){ invisible(unlist(as.list(x))) }
tstat <- function(r,n){ invisible(r * sqrt((n-2) / (1 - r*r))) } 
toZ <- function(r){ invisible(0.5 * log((1 + r) / (1 - r)))}
#toZ <- function(r){ invisible(atanh(r)) }

estimate <- function(val, permutations){
  invisible(1.0 - (getidx(val, permutations, length(permutations))/ length(permutations)))
}

getidx <- function(val, permutations, nperms){
  for(i in 1:(nperms-1)){ if(permutations[i] >= val) return(i); }
  return (nperms-1);
}

toP <- function(x,permutations){
  permutations <- sort(permutations)
  apply(x,c(1, 2),estimate,permutations)
}

plotInfo <- function(n1, n2, traits=100, n.perms=1000, strategy= c("Single","Slice","Box")){
  AA <- cor(matrix(rnorm(n1*traits),n1,traits))
  BB <- cor(matrix(rnorm(n2*traits),n2,traits))
  op <- par(mfrow=c(2,4))

  #Absolute difference in correlation
  plot(c(-0.5,0.5),c(-0.5,0.5),t='n',main="|CorA-CorB|")
  points(tolist(AA), tolist(BB), col=gray(1-abs(AA-BB)), pch=20, cex=0.1)
  abline(h=0,v=0,lty=3)
  
  #DCOR score
  DCOR <- (.5*(sign(AA)*AA^2-sign(BB)*BB^2))^2
  DCOR[is.nan(DCOR)] <- 0
  plot(c(-0.5,0.5),c(-0.5,0.5),t='n',main="DCOR scores (CorA-CorB)^2")
  points(tolist(AA), tolist(BB), col=gray(1-tolist(DCOR)/max(DCOR)), pch=20, cex=0.1)
  abline(h=0,v=0,lty=3)

  perms <- NULL
  while(length(perms) < n.perms){
    AAp   <- cor(matrix(rnorm(n1*traits),n1,traits))
    BBp   <- cor(matrix(rnorm(n2*traits),n2,traits))
    DCORp <- (AAp-BBp)^2
    DCORp <- (.5*(sign(AAp)*AAp^2-sign(BBp)*BBp^2))^2
    DCORp[is.nan(DCORp)] <- 0
    if(strategy[1]=="Single") perms <- c(perms, tolist(DCORp))
    if(strategy[1]=="Slice")  perms <- c(perms, apply(DCORp,1,max))
    if(strategy[1]=="Box")    perms <- c(perms, max(DCORp))
  }
  perms <- perms[sample(length(perms),n.perms)]
  dlikelihood <- toP(DCOR,perms)

  #DCOR likelihood
  plot(c(-0.5,0.5),c(-0.5,0.5),t='n',main="DCOR likelihood P[ (CorA-CorB)^2 ]")
  points(tolist(AA), tolist(BB), col=gray(tolist(dlikelihood)), pch=20, cex=0.1)
  abline(h=0,v=0,lty=3)
  abline(0, 1,lty=3)

  #LUDE
  plot(c(-0.5,0.5),c(-0.5,0.5),t='n',main="likelihood P[ (CorA-CorB) ] (Method: LUDE)")
#  LCOR <- pnorm((qnorm(dt(tstat(AA,n1),n1+n2))-qnorm(dt(tstat(BB,n2),n1+n2))))
#  LCOR[is.nan(LCOR)] <- 0.5
#  LCOR <- (0.5-abs(LCOR - 0.5))*2
  SE     <- sqrt((1/(n1-3)) + (1/(n2-3)))
  LCOR   <- pnorm( (toZ(AA) - toZ(BB)) / SE)
  LCOR[is.nan(LCOR)] <- 0.5
  LCOR <- (0.5-abs(LCOR - 0.5))*2
  points(tolist(AA), tolist(BB), col=gray(tolist(LCOR)), pch=20, cex=0.1)
  abline(h=0,v=0,lty=3)
  abline(0, 1,lty=3)

  hist(abs(AA-BB),main="Histogram of |CorA-CorB|")
  hist(tolist(DCOR),main="Histogram of (CorA-CorB)^2")
  hist(tolist(dlikelihood),main="Histogram of P[ (CorA-CorB)^2 ]")
  hist(tolist(LCOR),main="Histogram of P[ (CorA-CorB) ]")
  list(tolist(AA), tolist(BB), DCOR, toP(DCOR,perms), tolist(LCOR))
}

res <- plotInfo(100,400,50)

#res <- plotInfo(100,100,str="Slice")
#res <- plotInfo(100,100,str="Box")

#res <- plotInfo(100,300)
#res <- plotInfo(100,300,str="Slice")
#res <- plotInfo(100,300,str="Box")
