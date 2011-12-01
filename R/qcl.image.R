#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Nov, 2010
# 
# Image plot routines for QCL analysis
#

image.QCLscan <- function(x, QCLpermute, qcl.threshold = 0.35, against = c("markers","phenotypes"), do.grid=TRUE, grid.col = "black", verbose = FALSE, ...){
  colorrange <- c("white",gray.colors(10)[10:1])
  if(missing(QCLpermute)){
    if(verbose) cat("No permutations, using a cut-off leading to sub-optimal results\n")
    mymatrix <- QCLprofiles(x, qcl.threshold=qcl.threshold, against=against)
    mainlabel <- paste("QCLs at",qcl.threshold,"phenotypes vs",against[1])
  }else{
    mymatrix <- NULL
    for(p in 1:length(x)){
      if(verbose) cat("Processing:",p,"from QCL to LOD\n")
      mymatrix <- rbind(mymatrix,QCLtoLODvector(x, QCLpermute, p, against))
    }
    rownames(mymatrix) <- unlist(lapply(x,attr,"name"))
    mainlabel <- paste("QCL phenotypes vs",against[1])
  }
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix),
          main=mainlabel,
          yaxt="n",xaxt="n",ylab="", xlab="",col=colorrange,cex.main=0.7,breaks = c(1,2,3,4,5,6,7,8,9,10,11,100))
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix),las=2,cex.axis=0.5)
    axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.5)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd=1)
      abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
  }
  invisible(mymatrix)
}
