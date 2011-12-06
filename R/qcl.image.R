#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Nov, 2010
# 
# Image plot routines for QCL analysis
#

image.QCLobject <- function(x, against = c("markers","phenotypes"), do.grid=TRUE, grid.col = "black", verbose = FALSE, ...){
  colorrange <- c("white",gray.colors(10)[10:1])
  mymatrix <- NULL
  mynames <- NULL
  for(p in 1:length(x)){
    if(verbose) cat("Processing:",p,"from QCL to LOD\n")
    mymatrix <- rbind(mymatrix,QCLtoLODvector(x[[p]], against))
    mynames <- c(mynames,attr(x[[p]]$s,"name"))
  }
  rownames(mymatrix) <- mynames
  mainlabel <- paste("QCL phenotypes vs",against[1])
  if(!is.null(mymatrix)){ 
    image(1:ncol(mymatrix),1:nrow(mymatrix),t(mymatrix),
          main=mainlabel,
          yaxt="n",xaxt="n",ylab="", xlab="",col=colorrange,cex.main=0.7,breaks = c(0,1,2,3,4,5,6,7,8,9,10,100))
    axis(2,rownames(mymatrix),at=1:nrow(mymatrix),las=2,cex.axis=0.5)
    axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.5)
    if(do.grid){
      abline(h=seq(-0.5,nrow(mymatrix)+0.5,1),col=grid.col,lwd=1)
      abline(v=seq(-0.5,ncol(mymatrix)+0.5,1),col=grid.col,lwd=1)
    }
    box()
  }
  invisible(mymatrix)
}
