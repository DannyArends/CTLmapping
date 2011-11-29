#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Oct, 2011
# first written Nov, 2010
# 
# Plotting routines for QCL analysis
#

plot.QCLscan <- function(x, pheno.col = 1, qcl.threshold =0.3, do.legend=TRUE, ...){
  npheno <- length(x)
  if(pheno.col > npheno) stop("No such phenotype")
  pname <- attr(x[[pheno.col]],"name") 
  if(max(QCLprofiles(x,qcl.threshold)[pname,]) == 0) stop(paste("Threshold too high"))
  totpheno <- dim(x[[pheno.col]])[1]
  totmarkers <- dim(x[[pheno.col]])[2]
  y_range <- c(0,1.25*max(QCLprofiles(x,qcl.threshold)[pname,]))
  plot(c(0,totmarkers),y_range,type="n",main=paste("QCL mapping of",attr(x[[pheno.col]],"name")),ylab="# of significant QCL", xlab="Genetic marker")
  colorz <- NULL
  max_value <- max(abs(x[[pheno.col]]))
  for(t in seq(qcl.threshold,max_value,0.05)){
    points(QCLprofiles(x,t)[pname,], lwd=2,col=rgb((1/max_value)*t,0,0),type='l')
    colorz <- c(colorz,rgb((1/max_value)*t,0,0))
  }
  if(do.legend){
    legend("topleft",paste("Threshold =",seq(qcl.threshold,max_value,0.05)),lwd=2,col=colorz)
  }
}

plotAsLOD <- function(QTLscores, QCLscan, permutations, pheno.col = 1){
  npheno <- length(QCLscan)
  if(pheno.col > npheno) stop("No such phenotype")
  pname <- paste("Comparison QCL:QTL of",attr(QCLscan[[pheno.col]],"name"))
  QCLscores <- QCLtoLODvector(QCLscan, permutations, pheno.col=pheno.col)
  plot(c(0,length(QCLscores)),c(0,max(c(QCLscores,QTLscores))),type='n', main=pname, ylab="LOD",xlab="Marker")
  points(QCLscores,type='l',col="black",lwd=3)
  points(QTLscores,type='l',col="red",lwd=2,lty=3)
  legend("topleft",c("QCL","QTL"),col=c("black","red"),lty=c(1,3),lwd=c(3,2))
}

plotAsStackedHist <- function(qcl_result, qcl_perms, pheno.col=1, onlySignificant = TRUE, do.legend=TRUE, ...){
  summarized <- QCLtoLODvector(qcl_result, qcl_perms, pheno.col=pheno.col)
  plot(summarized, type='l',main=paste("Phenotype contribution to QCL of",attr(qcl_result[[pheno.col]],"name")),...)
  p <- rep(0,ncol(qcl_result[[pheno.col]]))
  i <- 1;
  QCLmatrix <- QCLtoLOD(qcl_result, qcl_perms, pheno.col, onlySignificant)
  if(!onlySignificant && ncol(QCLmatrix) > 15){
    cat("Warning disabled the legend please use onlySignificant = TRUE")
    do.legend=FALSE
  }
  mycolors <- rainbow(ncol(QCLmatrix))
  apply(QCLmatrix,2,
    function(d){
     for(idx in 1:length(d)){
        rect(idx-0.5,p[idx],idx+0.5,p[idx]+d[idx],col=mycolors[i])
      }
      p <<- p + d
      i <<- i + 1
    }
  )
  if(do.legend) legend("topleft",colnames(QCLmatrix),col=mycolors,lwd=1,cex=0.7)
  points(summarized,type='l',lwd=2)
  invisible(QCLmatrix)
}
