plotTraits <- function(genotypes, phenotypes, pheno.col = c(1, 2), marker = 1, doLog = FALSE){
  t1 <- phenotypes[, pheno.col[1]]
  t2 <- phenotypes[, pheno.col[2]]
  if(doLog){ t1    <- log10(t1); t2    <- log10(t2) }

  geno  <- genotypes[,marker]

  betas <- stdSlopeEstimate(t1,t2, geno)
  inter <- stdSlopeIntercept(t1,t2, geno, betas)

  t1name <- colnames(phenotypes)[pheno.col[1]]
  t2name <- colnames(phenotypes)[pheno.col[2]]

  plot(t1, t2, type = "n", main = "Trait-Trait scatterplot", xlab=t1name, ylab=t2name)
  max1 <- max(t1*1.25, na.rm=TRUE);  max2 <- max(t2*1.25, na.rm=TRUE)
  min1 <- min(t1, na.rm=TRUE);       min2 <- min(t2, na.rm=TRUE)
  abline( v = seq(min1, max1, (max1 - min1)/20 ), lty = 3, col = colors()[ 440 ] )
  abline( h = seq(min2, max2, (max2 - min2)/20 ), lty = 3, col = colors()[ 440 ] )

  points(t1, t2, pch=19, cex=0.5, col = geno)
  
  cors <- NULL

  for(x in 1:length(unique(geno))){
    idx <- which(geno == unique(geno)[x])
    cors <- c(cors, cor(t1[idx], t2[idx]))
    abline(inter[x], betas[x], col = x, lwd=1.5)
  }
  legend("topright", title="Correlation", col=unique(geno), legend = round(cors, d=2),lwd=1, cex=0.7)
  legend("topleft",  title="Slope", col=unique(geno), legend = round(cors, d=2),lwd=1, cex=0.7)
}



x <- seq(1,24)


p <- function(genotypes, phenotypes, pheno.col=c(1,2), doLOG=FALSE){
  nmar <- ncol(genotypes)
  t1coords <- as.numeric(unlist(phenotypes[, pheno.col[1]]))
  t1coords <- as.numeric(unlist(phenotypes[, pheno.col[2]]))
  col <- as.numeric(unlist(genotypes))
  x <- rep(t1coords, nmar)
  y <- unlist(lapply(seq(1,nmar), rep, length(t1coords)))
  z <- rep(t2coords, nmar)
  if(doLOG) x <- log10(x)
  if(doLOG) y <- log10(y)
  scatterplot3d(x,y,z, color=col, pch=0.2)
}


p <- function(ctls, doLOG=FALSE){
  nmar <- ncol(genotypes)
  z <- unlist(lapply(ctls,function(x){
    as.numeric(unlist(x$ctl))
  }))
  id <- which(z == 0)
  nt <- (length(z)/length(ctls))/nmar
  x <- unlist(lapply(seq(1,nt), rep, nmar))
  x <- rep(x,nt)
  y <- unlist(lapply(seq(1,nt), rep, nt))
  y <- rep(y,nmar)

  col <- unlist(lapply(seq(1,nmar), rep, nt))
  col <- rep(col, nt)

  scatterplot3d(x[-id],y[-id],z[-id], color=col[-id], pch=19, cex.symbols=0.1)

  t1coords <- as.numeric(unlist(phenotypes[, pheno.col[1]]))
  t1coords <- as.numeric(unlist(phenotypes[, pheno.col[2]]))
#  col <- as.numeric(unlist(genotypes))
  x <- rep(t1coords, nmar)
  y <- unlist(lapply(seq(1,nmar), rep, length(t1coords)))
  z <- rep(t2coords, nmar)
  if(doLOG) x <- log10(x)
  if(doLOG) z <- log10(z)
  scatterplot3d(x,y,z, pch=19, cex.symbols=0.1)
}

test <-
n <- 300

qtl1 <- 3 ; noise1 <- 2
qtl2 <- 0 ; noise2 <- .001

t1aa <-  .5*qtl1 + (rnorm(n) * noise1)
t1ab <- sort(rnorm(2*n)) * noise1*rnorm(2*n)
t1bb <- -.5*qtl1 + (rnorm(n) * noise1)

t2aa <-  .5*qtl2 + (rnorm(n) * noise2)
t2ab <- sort(rnorm(2*n)) * noise2*rnorm(2*n)
t2bb <- -.5*qtl2 + (rnorm(n) * noise2)

t1 <- c(t1aa, t1ab, t1bb)
t2 <- c(t2aa, t2ab, t2bb)

geno <- c(rep(1,n), rep(2,2*n), rep(3,n))

p(t(t(geno)), cbind(t1,t2))

