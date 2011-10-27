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
  if(max(QCLscanToProfile(x[[pheno.col]],qcl.threshold)) == 0) stop(paste("Threshold too high"))
  totpheno <- dim(x[[pheno.col]])[1]
  totmarkers <- dim(x[[pheno.col]])[2]
  y_range <- c(0,1.25*max(QCLscanToProfile(x[[pheno.col]],qcl.threshold)))
  plot(c(0,totmarkers),y_range,type="n",main=paste("QCL mapping of ",attr(x[[pheno.col]],"name")),ylab="# of significant QCL", xlab="Genetic marker")
  colorz <- NULL
  for(t in seq(qcl.threshold,max(x[[pheno.col]]),0.05)){
    points(QCLscanToProfile(x[[pheno.col]],t), lwd=2,col=rgb((1/max(x[[pheno.col]]))*t,0,0),type='l')
    colorz <- c(colorz,rgb((1/max(x[[pheno.col]]))*t,0,0))
  }
  legend("topleft",paste("Threshold =",seq(qcl.threshold,max(x[[pheno.col]]),0.05)),lwd=2,col=colorz)
}

image.QCLscan <- function(x, qcl.threshold =0.35, against = c("markers","phenotypes"), ...){
  colorrange <- c("white",gray.colors(100)[100:1])
  if(against[1] == "markers"){
    mymatrix <- NULL
    ylabels <- NULL
    for(p in 1:length(x)){
      mymatrix <- rbind(mymatrix,QCLscanToProfile(x[[p]],qcl.threshold))
      ylabels <- c(ylabels,attr(x[[p]],"name"))
    }
    image(1:ncol(mymatrix),1:length(x),t(mymatrix),main=paste("QCLs at",qcl.threshold),yaxt="n",yaxt="n",ylab="", xlab="Genetic marker",col=colorrange)
    axis(2,ylabels,at=1:length(x),las=2,cex.axis=0.7)
    invisible(return(mymatrix))
  }
  if(against[1] == "phenotypes"){
    phenotypenames <- NULL
    ylabels <- NULL
    for(p in 1:length(x)){
      phenotypenames <- unique(c(phenotypenames,unique(as.character(unlist(apply(x[[p]],2,function(x){names(which(x > qcl.threshold))}))))))
      ylabels <- c(ylabels,attr(x[[p]],"name"))
    }
    mymatrix <- matrix(0,length(x),length(phenotypenames))
    colnames(mymatrix) <- phenotypenames
    for(p in 1:length(x)){
      current_table <- table(as.character(unlist(apply(x[[p]],2,function(x){names(which(x > qcl.threshold))}))))
      mymatrix[p,names(current_table)] <- current_table
    }
    image(1:ncol(mymatrix),1:length(x),t(mymatrix),main=paste("QCLs at",qcl.threshold),yaxt="n",xaxt="n",ylab="", xlab="",col=colorrange)
    axis(1,colnames(mymatrix),at=1:ncol(mymatrix),las=2,cex.axis=0.7)
    axis(2,ylabels,at=1:length(x),las=2,cex.axis=0.7)
    invisible(return(mymatrix))
  }
}
