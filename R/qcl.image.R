#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Nov, 2010
# 
# Image plot routines for QCL analysis
#

image.QCLscan <- function(x, qcl.threshold = 0.35, against = c("markers","phenotypes"), do.grid=TRUE, grid.col = "black", ...){
  colorrange <- c("white",gray.colors(100)[100:1])
  mymatrix <- QCLprofiles(x, qcl.threshold=qcl.threshold, against=against)
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix),
          main=paste("QCLs at",qcl.threshold,"phenotypes vs",against[1]),
          yaxt="n",xaxt="n",ylab="", xlab="",col=colorrange,cex.main=0.7)
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix),las=2,cex.axis=0.5)
    axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.3)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd=1)
      abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
  }
  invisible(mymatrix)
}
