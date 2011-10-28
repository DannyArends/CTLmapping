#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Nov, 2010
# 
# Plotting routines for QCL analysis
#

plot.QCLscan <- function(x, pheno.col = 1, qcl.threshold =0.3, ...){
  npheno <- length(x)
  if(pheno.col > npheno) stop("No such phenotype")
  pname <- attr(x[[pheno.col]],"name") 
  if(max(QCLprofiles(x,qcl.threshold)[pname,]) == 0) stop(paste("Threshold too high"))
  totpheno <- dim(x[[pheno.col]])[1]
  totmarkers <- dim(x[[pheno.col]])[2]
  y_range <- c(0,1.25*max(QCLprofiles(x,qcl.threshold)[pname,]))
  plot(c(0,totmarkers),y_range,type="n",main=paste("QCL mapping of ",attr(x[[pheno.col]],"name")),ylab="# of significant QCL", xlab="Genetic marker")
  colorz <- NULL
  for(t in seq(qcl.threshold,max(x[[pheno.col]]),0.05)){
    points(QCLprofiles(x,t)[pname,], lwd=2,col=rgb((1/max(x[[pheno.col]]))*t,0,0),type='l')
    colorz <- c(colorz,rgb((1/max(x[[pheno.col]]))*t,0,0))
  }
  legend("topleft",paste("Threshold =",seq(qcl.threshold,max(x[[pheno.col]]),0.05)),lwd=2,col=colorz)
}
