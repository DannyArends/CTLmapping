#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Nov, 2010
# 
# Plotting routines for QCL analysis
#

print.QCLscan <- function(x, ...){
  cat("QCLscan summary\n\n")
  cat("- Number of phenotypes scanned",length(x),"/",dim(x[[1]])[1],"\n")
  cat("- Number of markers",dim(x[[1]])[2],"\n")
  unlist(...)
}

plot.QCLscan <- function(x, pheno.col = 1, qcl.threshold =0.3, ...){
  npheno <- length(x)
  if(pheno.col > npheno) stop("No such phenotype")
  if(max(QCLprofiles(x[[pheno.col]],qcl.threshold)) == 0) stop(paste("Threshold too high"))
  totpheno <- dim(x[[pheno.col]])[1]
  totmarkers <- dim(x[[pheno.col]])[2]
  y_range <- c(0,1.25*max(QCLprofiles(x[[pheno.col]],qcl.threshold)))
  plot(c(0,totmarkers),y_range,type="n",main=paste("QCL mapping of ",attr(x[[pheno.col]],"name")),ylab="# of significant QCL", xlab="Genetic marker")
  colorz <- NULL
  for(t in seq(qcl.threshold,max(x[[pheno.col]]),0.05)){
    points(QCLprofiles(x[[pheno.col]],t), lwd=2,col=rgb((1/max(x[[pheno.col]]))*t,0,0),type='l')
    colorz <- c(colorz,rgb((1/max(x[[pheno.col]]))*t,0,0))
  }
  legend("topleft",paste("Threshold =",seq(qcl.threshold,max(x[[pheno.col]]),0.05)),lwd=2,col=colorz)
}

image.QCLscan <- function(x, qcl.threshold = 0.35, against = c("markers","phenotypes"), do.grid=TRUE, grid.col = "black", ...){
  colorrange <- c("white",gray.colors(100)[100:1])
  mymatrix <- QCLprofiles(x, qcl.threshold=qcl.threshold, against=against)
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix),
          main=paste("QCLs at",qcl.threshold,"phenotypes vs",against[1]),
          yaxt="n",xaxt="n",ylab="", xlab="",col=colorrange)
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix),las=2,cex.axis=0.7)
    axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.7)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd=1)
      abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
  }
  invisible(mymatrix)
}
