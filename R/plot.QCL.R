#
#
# plot.QCL.R
#
# copyright (c) 2010 Danny Arends and Bruno Tesson
# last modified Jun, 2011
# first written nov, 2010
# 
# Plotting routines for QCL analysis
#

#Heatmap the output of a markerQCL object
#Returns the hclust object used to order the traits shown in the heatmap
plot.markerQCL <- function(markerQCL, QCL.threshold=0.5, QCL.significant = 0, ...){
  aboveThreshold <- count.QCL(markerQCL, QCL.threshold)
  difCorrelated <- which(aboveThreshold > QCL.significant)
  if(length(difCorrelated) <= 1){
     warning("No phenotype shows differential correlation with more than: ",QCL.significant," other phenotypes at QCL.threshold: ",QCL.threshold,"\n")
     return()
  }
  ccorclass1 <- markerQCL[[2]][difCorrelated,difCorrelated]
  ccorclass2 <- markerQCL[[3]][difCorrelated,difCorrelated]
  if(nrow(ccorclass1) >= 2){
    ordering <- hclust(dist(ccorclass1))
    ccorclass1 <- ccorclass1[ordering$order,ordering$order]
    ccorclass2 <- ccorclass2[ordering$order,ordering$order]
    upper <- t(upper.tri(ccorclass1))*ccorclass1
    lower <- t(lower.tri(ccorclass2))*ccorclass2
    op <- par(mar=c(8,8,4,2)+0.1)
    genotypes <- upper.tri(matrix(0,nrow(ccorclass1),nrow(ccorclass1)))
    for(x in 1:nrow(genotypes)){
      genotypes[x,x] <- NA
    }
    image(seq(0,nrow(ccorclass1)),seq(0,nrow(ccorclass1)),genotypes,ylab="",xlab="",yaxt="n",xaxt="n",col=c(rgb(255,0,0,50,maxColorValue=255),rgb(255,255,0,50,maxColorValue=255)),main=paste("QCL at marker",attr(markerQCL,"marker")),...)
    image(seq(0,nrow(ccorclass1)),seq(0,nrow(ccorclass1)),lower+upper,breaks=c(-1,-0.75,-0.5,0.5,0.75,1),col=c("blue","lightblue",rgb(255,255,255,155,maxColorValue=255),"lightgreen","green"),add=TRUE,...)
    axis(2,at=seq(0.5,nrow(ccorclass1)),labels=colnames(ccorclass1),las=1)
    axis(1,at=seq(0.5,nrow(ccorclass1)),labels=colnames(ccorclass1),las=3)
    if(nrow(genotypes) < 50){
      grid(nrow(ccorclass1),nrow(ccorclass1),lwd=1,lty=1,col="black")
    }else{
      cat("INFO: ",attr(markerQCL,"marker"),"- 50+ traits, plot grid disabled, Add with: grid(",nrow(genotypes),",",nrow(genotypes),",lwd=1,lty=1,col=\"black\")\n")
    }
    abline(v=0)
    abline(h=nrow(ccorclass1))
    ordering
  }else{
    warning("No phenotype shows differential correlation with more than: ",QCL.significant," other phenotypes at QCL.threshold: ",QCL.threshold,"\n")
  }
}

#Heatmap the output of a markerQCL object (almost the same as above)
plot.markerQCL2 <- function(markerQCL, peekheight=1){
  selection <- which(markerQCL[[4]] >= peekheight)
  if(length(selection) > 0){
    g1 <- markerQCL[[2]][selection,selection]
    g2 <- markerQCL[[3]][selection,selection]
    ordering <- markerQCL[[1]][selection,selection]
    clustering <- hclust(dist(ordering))

    upper <- upper.tri(g1)*sign(g1)*(g1[clustering$order,clustering$order]^2)
    lower <- lower.tri(g2)*sign(g2)*(g2[clustering$order,clustering$order]^2)
    colorz <- c("red","yellow","white","lightblue","blue")
    heatmap(t(upper+lower),col=colorz,breaks=c(-1,-0.5,-0.3,0.3,0.5,1),Colv=NA,Rowv=NA,scale="none",main=paste("Correlation at: ",attr(markerQCL,"marker")))
    clustering
  }else{
    plot(1:10)
    return(NULL)
  }
}

#Heatmap the output of a QCLscan object
image.QCLscan <- function(QCLscan, QCL.thresholds=0, cluster=FALSE){
  if(cluster){
    heatmap(t(QCLscan[,which(apply(QCLscan,2,function(x){max(x) > QCL.thresholds}))]),Colv=NA,scale="none",col=c("white",gray.colors(10)[10:1],"black"))
  }else{
    heatmap(t(QCLscan[,which(apply(QCLscan,2,function(x){max(x) > QCL.thresholds}))]),Colv=NA, Rowv=NA,scale="none",col=c("white",gray.colors(10)[10:1],"black"))
  }
}

getMarkerOffsets <- function(cross, cmBetween=25){
  offsets <- unlist(lapply(pull.map(cross),max))
  offsets <- offsets+cmBetween
  offsets <-c(0,offsets)

  cnt <- 1
  myoffsets <- NULL
  for(x in nmar(cross)){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(unlist(pull.map(cross)))
  mlocations
}

bplot.QCLscan <- function(cross, QCLscan,addQTL=TRUE){
  phenoname <- attr(QCLscan,"phenoname")
  offsets <- getMarkerOffsets(cross)
  plot(offsets,as.numeric(c(max(QCLscan),rep(min(QCLscan),sum(nmar(cross))-1))),type="n",main=phenoname,xlab="",ylab="")
  
  olnm <- 1
  snm <- 0
  if(addQTL){
    ids <- which(phenames(cross) == phenoname)
    cat(phenoname,"=",ids,"\n")
    cross <- calc.genoprob(cross)
    qtlscan <- scanone(cross, pheno.col=ids)
    dtmax <- max(QCLscan)
    qtlmax <- max(qtlscan[,3],10)
    qtlscan[,3] <- qtlscan[,3]*(dtmax/qtlmax)
    cat(qtlmax," ",dtmax,"\n")
    axis(4,at=seq(0,1.7*dtmax,1),round((qtlmax/dtmax) * seq(0,1.7*dtmax,1),1))
  }
  for(nm in nmar(cross)){
    snm = nm + snm
    points(offsets[olnm:snm],as.numeric(qtlscan[olnm:snm,3]),lwd=3,type='o',col="red")
    points(offsets[olnm:snm],as.numeric(QCLscan[olnm:snm,phenoname]),lwd=2,type='o',col="blue")
    olnm=snm+1
  }
  boxplot(t(QCLscan),at=offsets,main=phenoname,add=T,xlab=F,names=rep("",length(offsets)),ylab="QCL",pch=20)
  legend("topright",c("Boxplot QCL","QCL","QTL (scanone)"),lty=c(1,1,1),lwd=c(1,3,2),col=c("black","blue","red"))
}


#plot.QCLscan Plot a QCL profile of a single phenotype
#Optionally add a scanone object to overlay the QTL profile
plot.QCLscan <- function(cross, QCLscan, QCL.thresholds=NULL, addQTL=FALSE, ...){
  todrop <- markernames(cross)[unlist(1:sum(nmar(cross)))[-attr(QCLscan,"markers")]]
  cross <- drop.markers(cross,todrop)
  phenoname <- attr(QCLscan,"phenoname")
  QCLscanProfile <- lodscorestoscanone(cross,apply(QCLscan,1,var),traitnames = "QCL")
  if(addQTL){
    ids <- which(phenames(cross) == phenoname)
    cat(phenoname,"=",ids,"\n")
    cross <- calc.genoprob(cross)
    qtlscan <- scanone(cross, pheno.col=ids)
  }
  dtmax <- max(QCLscanProfile[,3])
  if(addQTL) qtlmax <- max(qtlscan[,3],10)
  if(addQTL) QCLscanProfile[,3] <- QCLscanProfile[,3]*(qtlmax/dtmax)
  if(addQTL){
    op <- par(mar=c(5, 4, 4, 5) + 0.1)
    plot(qtlscan,QCLscanProfile,y=c(0,1.7*qtlmax),col=c("red","black"),lty=c(1,1),lwd=c(3,2),main=phenoname,...)
    axis(4,at=seq(0,1.7*qtlmax,1),round((dtmax/qtlmax) * seq(0,1.7*qtlmax,1),1))
    legend("topright",c("QTL (scanone)","QCL"),lty=c(1,1),lwd=c(3,2),col=c("red","black"))
  }else{
    plot(QCLscanProfile,y=c(0,1.7*dtmax),col="black",main=phenoname,ylab="QCL",...)
    if(!is.null(QCL.thresholds)){
      colorz <- c("red","orange","green")
      i <- 1
      for(x in QCL.thresholds){
        abline(h=x,col=colorz[i],lty=3,lwd=2)
        i <- i+1
      }
      legend("topright",c("QCL",names(QCL.thresholds)),lwd=2,col=c("black",colorz),lty=c(1,3,3,3))
    }else{
      legend("topright",c("QCL"),lwd=1,col=c("black"))
    }
  }
  if(addQTL) QCLscanProfile[,3] <- QCLscanProfile[,3]*(dtmax/qtlmax)
  invisible(QCLscanProfile)
}

#Plot the detailed correlations and QCL score of a given phenotype in a markerQCL object, 
#when threshold is specified only the the elements above the threshold are plotted
plot.markerQCL.detailed <- function(markerQCL, pheno.col=1, difCorThreshold=0){
  sorted <- sort(abs(markerQCL[[1]][,pheno.col]),decreasing=TRUE)
  sorted <- sorted[which(sorted > difCorThreshold)]
  ordering <- names(sorted)
  plot(c(1,length(ordering)),c(0,1),xlab="Other phenotypes",ylab="Correlation",main=paste("Correlation probe ",pheno.col," at marker",attr(markerQCL,"marker"),sep=" "),type='n')
  points(markerQCL[[2]][ordering,pheno.col]^2,pch=20,col='red',type='s')
  points(markerQCL[[3]][ordering,pheno.col]^2,pch=20,col='green',type='s')
  points(abs(markerQCL[[1]][ordering,pheno.col]),pch=20,col='blue',type='l',lwd=4)
  legend("topright",c("Genotype 1","Genotype 2","Absolute difference"),pch=20,col=c("red","green","blue"))
  invisible(ordering)
}
