### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: article.Rnw:90-91
###################################################
options(width=87, digits=3, scipen=4)


###################################################
### code chunk number 2: article.Rnw:94-97
###################################################
  require(ctl)
  set.seed(200)
  data(ath.result)


###################################################
### code chunk number 3: permutation_scores
###################################################
op <- par(cex=0.5)
hist(ath.result,c(3,2,1))
op <- par(cex=1)


###################################################
### code chunk number 4: singletrait
###################################################
op <- par(cex=0.5)
plot(ath.result[[1]])
op <- par(cex=1)


###################################################
### code chunk number 5: singletraitLOD
###################################################
op <- par(cex=0.5)
CTLasLOD(ath.result[[1]], ath.result[[1]]$qtl)
op <- par(cex=1)


###################################################
### code chunk number 6: p2mmatrix
###################################################
par(mar=c(3,8,1.1,1.1))
plot(ath.result, against="markers", do.grid=TRUE, grid.col="white")


###################################################
### code chunk number 7: QTLmatrix
###################################################
par(mar=c(3,8,1.1,1.1))
QTLimage(ath.result, do.grid=TRUE, grid.col="white")


###################################################
### code chunk number 8: p2pmatrix
###################################################
par(mar=c(3,8,1.1,1.1))
plot(ath.result, against="phenotypes", do.grid=TRUE, grid.col="white")


