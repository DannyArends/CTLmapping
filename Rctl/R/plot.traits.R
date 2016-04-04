#
# plot.traits.R
#
# copyright (c) 2010-2012 - GBIC, Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified Apr, 2013
# first written Apr, 2013
# 
# Trait vs Trait scatter plot routine for CTL analysis
#

plotTraits <- function(genotypes, phenotypes, phenocol = c(1, 2), marker = 1, doRank = FALSE){
  t1 <- phenotypes[, phenocol[1]]
  t2 <- phenotypes[, phenocol[2]]
  if(doRank){ t1 <- rank(t1); t2 <- rank(t2) }

  geno  <- genotypes[,marker]

  betas <- stdSlopeEstimate(t1,t2, geno)
  inter <- stdSlopeIntercept(t1,t2, geno, betas)

  t1name <- colnames(phenotypes)[phenocol[1]]
  t2name <- colnames(phenotypes)[phenocol[2]]

  plot(t1, t2, type = "n", main = "Trait-Trait scatterplot", xlab=t1name, ylab=t2name)
  max1 <- max(t1*1.25, na.rm=TRUE);  max2 <- max(t2*1.25, na.rm=TRUE)
  min1 <- min(t1, na.rm=TRUE);       min2 <- min(t2, na.rm=TRUE)
  abline( v = seq(min1, max1, (max1 - min1) / 40), lty = 3, col = colors()[ 440 ] )
  abline( h = seq(min2, max2, (max2 - min2) / 40), lty = 3, col = colors()[ 440 ] )

  points(t1, t2, pch=19, cex=0.5, col = geno)
  cors <- NULL
  for(x in 1:length(unique(geno))){
    idx <- which(geno == unique(geno)[x])
    cors <- c(cors, cor(t1[idx], t2[idx],use="pair"))
    abline(inter[x], betas[x], col = x, lwd=1.5)
  }
  legend("topright", title="Correlation", col=unique(geno), legend = round(cors, digits = 2), lwd=1, cex=0.7)
  legend("topleft",  title="Slope", col=unique(geno), legend = round(betas, digits = 2), lwd=1, cex=0.7)
}

# end of plot.traits.R

