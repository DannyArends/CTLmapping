#
# qcl.plot.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jan, 2012
# first written Nov, 2010
# 
# Plotting routines for QCL analysis
#

QCLasLOD <- function(QCLscan, QTLscores, main, do.legend=TRUE){
  if(missing(QCLscan)) stop("argument 'QCLscan' is missing, with no default")
  if(missing(main)){
    main <- paste("Comparison QCL:QTL of",attr(QCLscan$qcl,"name"))
  }
  op <- par(mfrow=c(2,1))
  plot.QCLscan(QCLscan,do.legend=do.legend)
  QCLscores <- QCLtoLODvector(QCLscan)
  if(!missing(QTLscores)){
    QTLscores <- as.numeric(unlist(QTLscores))
    plot(c(0,length(QCLscores)),c(0,max(c(QCLscores,QTLscores))),type='n', main=main, ylab="LOD",xlab="Marker")
    points(QCLscores,type='l',col="black",lwd=3)
    points(QTLscores,type='l',col="red",lwd=2,lty=1)
    if(do.legend) legend("topleft",c("QCL","QTL"),col=c("black","red"),lty=c(1,1),lwd=c(3,2))
  }else{
    plot(c(0,length(QCLscores)),c(0,max(c(QCLscores))),type='n', main=main, ylab="LOD",xlab="Marker")
    points(QCLscores,type='l',col="black",lwd=3)
  }
}

plot.QCLobject <- function(x, pheno.col=1:length(x), ...){
  if(length(x) == 1){
    return(plot.QCLscan(x[[1]],...))
  }
  if(length(pheno.col) == 1){
    return(plot.QCLscan(x[[pheno.col]],...))
  }
  return(image.QCLobject(x,...))
}

plot.QCLscan2 <- function(x, addQTL = TRUE, onlySignificant = TRUE, significance = 0.05, do.legend=TRUE, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  if(!is.null(x$p) && !is.nan(getPermuteThresholds(x,significance)[1])){
    mysign <- as.numeric(which(apply(abs(x$qcl),1,max) > getPermuteThresholds(x,significance)[1]))
  }else{
    mysign <- as.numeric(which(apply(abs(x$l),1,max) > -log10(significance)))
  }
  if(length(mysign) ==0){
    mysign <- 1:nrow(x$qcl)
    do.legend=FALSE
  }
  if(is.null(x$l)){
    QCLmatrix <- matrix(abs(x$qcl[mysign, ]),length(mysign),ncol(x$qcl))
  }else{
    QCLmatrix <- matrix(x$l[mysign, ],length(mysign),ncol(x$l))
  }
  summarized <- apply(QCLmatrix,2,sum)
  plot(c(0,ncol(x$qcl)),c(0,max(c(summarized,x$qtl))), type='n',xlab="Marker", ylab="-log10(P-value)", main=paste("Phenotype contribution to QCL of",attr(x$qcl,"name")),...)
  p <- rep(0,ncol(x$qcl))
  i <- 1;
  mycolors <- rep("black",nrow(QCLmatrix))
  apply(QCLmatrix,1,
    function(d){
      points(d,type='b',col=mycolors[i],pch=i,lwd=2)
      p <<- p + d
      i <<- i + 1
    }
  )
  if(do.legend){
    legend("topleft",c("QTL",paste("CTL",rownames(x$qcl)[mysign])),col=c("red",mycolors),lwd=2,pch=c(NA,1:length(mycolors)),cex=1.2)
  }
  points(as.numeric(x$qtl),type='l',lwd=2,col="red")
  rownames(QCLmatrix) <- rownames(x$qcl)[mysign]
  invisible(QCLmatrix)
}

plot.QCLscan <- function(x, addQTL = TRUE, onlySignificant = TRUE, significance = 0.05, do.legend=TRUE, ...){
  if(missing(x)) stop("argument 'x' is missing, with no default")
  if(!is.null(x$p) && !is.nan(getPermuteThresholds(x,significance)[1])){
    mysign <- as.numeric(which(apply(abs(x$qcl),1,max) > getPermuteThresholds(x,significance)[1]))
  }else{
    mysign <- as.numeric(which(apply(abs(x$l),1,max) > -log10(significance)))
  }
  if(length(mysign) ==0){
    mysign <- 1:nrow(x$qcl)
    do.legend=FALSE
  }
  if(is.null(x$l)){
    QCLmatrix <- matrix(abs(x$qcl[mysign, ]),length(mysign),ncol(x$qcl))
  }else{
    QCLmatrix <- matrix(x$l[mysign, ],length(mysign),ncol(x$l))
  }
  summarized <- apply(QCLmatrix,2,sum)
  plot(c(0,ncol(x$qcl)),c(0,max(c(summarized,x$qtl))), type='n',xlab="Marker", ylab="-log10(P-value)", main=paste("Phenotype contribution to QCL of",attr(x$qcl,"name")),...)
  p <- rep(0,ncol(x$qcl))
  i <- 1;
  mycolors <- topo.colors(nrow(QCLmatrix))
  apply(QCLmatrix,1,
    function(d){
     for(idx in 1:length(d)){
        rect(idx-0.5,p[idx],idx+0.5,p[idx]+d[idx],col=mycolors[i],lwd=0,lty=0)
      }
      p <<- p + d
      i <<- i + 1
    }
  )
  if(!is.null(x$l)){
    n <- dim(x$qcl)[2]
    abline(h=-log10(c(0.05/n,0.01/n,0.001/n)),col=c("red","orange","green"),lty=2)
    legend("topright",as.character(paste("QCL-FDR:",c(0.05,0.01,0.001),"%")),col=c("red","orange","green"),lty=rep(2,3),lwd=1,cex=0.7)
  }
  if(do.legend){
    legend("topleft",rownames(x$qcl)[mysign],col=mycolors,lwd=1,cex=0.7)
  }
  points(summarized,type='l',lwd=1)
  points(as.numeric(x$qtl),type='l',lwd=2,col="red")
  rownames(QCLmatrix) <- rownames(x$qcl)[mysign]
  invisible(QCLmatrix)
}

plot.QCLpermute <- function(x, type="s", ...){
  if(missing(x)) stop("argument 'x' which expects a 'QCLpermute' object is missing, with no default")
  plot(seq(0,0.9,0.01),QCLscoretoPvalue(seq(0,0.9,0.01),x),main="QCL to P.value",xlab="QCL",ylab="Pvalue", type=type,...)
  significant <- print.QCLpermute(x)
  mycolors <- c("red","orange","green")
  idx <- 1
  for(y in c(.05,.01,.001)){
    lines(rbind(c(-1,y),c(significant[idx],y)),lty=2,col=mycolors[idx])
    lines(rbind(c(significant[idx],-1),c(significant[idx],y)),lty=2,col=mycolors[idx])
    idx <- idx+1
  }
  legend("topright",c("QCL-FDR: 5%","QCL-FDR: 1%","QCL-FDR: 0.1%"), col=mycolors, lty=2, lwd=1, cex=0.7)
}

# end of qcl.plot.R
