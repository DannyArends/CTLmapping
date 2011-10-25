#QCL internal routine - Optimized using snow and 2 cores
#Does a Single marker return a list with:
# [[1]] difCorMatrix 
# [[2]] Correlation matrix individuals with genotype 1
# [[3]] Correlation matrix individuals with genotype 2
# [[4]] A vector of counted scores for each phenotype above the difCorThreshold
#Note: Also saves the object to: output/difCor<marker>.Rdata
QCLscan.internal.deprecated <- function(cross, marker, QCL.threshold=0.25, method="pearson", directory="output", saveRdata=FALSE, cpu_cluster){
  require(snow)
  expressions <- matrix(unlist(pull.pheno(cross)),dim(pull.pheno(cross))[1],dim(pull.pheno(cross))[2])
  colnames(expressions) <- colnames(pull.pheno(cross))
  genotypes <- pull.geno(cross)
  markerName <- markernames(cross)[marker]
  
  
  work <- vector("list",2)
  work[[1]] <- which(genotypes[,marker]==1)
  work[[2]] <- which(genotypes[,marker]==2)
  
  if(!is.null(cpu_cluster)){
    results <- parLapply(cpu_cluster,work,correlationmatrix,expressions=expressions,method=method)
  }else{
    results <- lapply(work,correlationmatrix,expressions=expressions,method=method)
  }
  
  QCL <- vector("list",4)
  QCL[[1]] <- 0.5*((sign(results[[1]]) * results[[1]]^2) -  (sign(results[[2]]) * results[[2]]^2))
  QCL[[2]] <- results[[1]]
  QCL[[3]] <- results[[2]]
  QCL[[4]] <- count.QCL(QCL,QCL.threshold)
  
  traitnames <- colnames(expressions)
  colnames(QCL[[1]]) <- traitnames
  rownames(QCL[[1]]) <- traitnames
  colnames(QCL[[2]]) <- traitnames
  rownames(QCL[[2]]) <- traitnames
  colnames(QCL[[3]]) <- traitnames
  rownames(QCL[[3]]) <- traitnames
  
  attr(QCL,"marker") <- markerName
  attr(QCL,"phenotypes") <- traitnames
  
  if(saveRdata) save(QCL,file=paste(directory,"/QCL_",marker,".Rdata",sep=""))
  class(QCL) <- c(class(QCL),"markerQCL")
  QCL
}

#Main routine to do the entire analysis
#Note: Does all the markers one by one (optimized to use 2 cores)
#Note: The difCor object in memory is very large
#Note: Based on the amount of traits and markers this could take a LONG time
QCLscan.deprecated <- function(cross, pheno.col=1, marker.col, QCL.threshold=0.25, significant = 0, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, saveRdata=FALSE, snow=TRUE, verbose=TRUE){
  s <- proc.time()
  if(doplot && !file.exists(directory)) dir.create(directory)
  if(verbose) cat("Analysis of ",ncol(cross$pheno)," traits at ",sum(nmar(cross))," markers\n")
  phenoname <- phenames(cross)[pheno.col]
  cross <- getCorrelatedPhenotypes(cross,pheno.col=pheno.col)
  if(nphe(cross)==0) stop("This phenotype doesn't show covariation with an other phenotype")
  if(!missing(marker.col)){
    totmarkers <- marker.col
  }else{
    totmarkers <- 1:sum(nmar(cross))
  }
  QCLscan <- NULL
  if(snow){
    require(snow)
    cpu_cluster <- makeCluster(c("localhost","localhost"))
  }else{
    cpu_cluster <- NULL
  }
  for(marker in totmarkers){
    sl <- proc.time()
    if(snow){
      gcLoop();
    }
    results <- QCLscan.internal.deprecated(cross, marker, QCL.threshold, method,directory,saveRdata, cpu_cluster)
    if(snow){
      gcLoop();
    }
    if(doplot){
      png(paste(directory,"/",marker,".jpg",sep=""))
        plot.markerQCL(results, QCL.threshold, significant)
      dev.off()
    }
    
    QCLscan <- rbind(QCLscan,results[[4]])
    
    results <- NULL
    if(verbose && marker %% 25 == 0){
      el <- proc.time()
      cat("Marker ",marker,"/ [",min(totmarkers),"..",max(totmarkers),"] took: ",as.numeric(el[3]-sl[3]),", so far:",as.numeric(el[3]-s[3]),"Seconds.\n")
    }
  }
  if(snow) stopCluster(cpu_cluster)
  cat("Analysis took: ",as.numeric(el[3]-s[3]),"Seconds\n")
  rownames(QCLscan) <- markernames(cross)[totmarkers]
  colnames(QCLscan) <- phenames(cross)
  class(QCLscan) <- c(class(QCLscan),"QCLscan")
  attr(QCLscan,"markers") <- totmarkers
  attr(QCLscan,"phenoname") <- phenoname
  if(writefile) write.table(QCLscan,file="eQCL_results.txt",sep="\t")
  QCLscan
}


#Heatmap the output of a markerQCL object
#Returns the hclust object used to order the traits shown in the heatmap
plot.markerQCL <- function(x, QCL.threshold=0.5, QCL.significant = 0, ...){
  aboveThreshold <- count.QCL(x, QCL.threshold)
  difCorrelated <- which(aboveThreshold > QCL.significant)
  if(length(difCorrelated) <= 1){
     warning("No phenotype shows differential correlation with more than: ",QCL.significant," other phenotypes at QCL.threshold: ",QCL.threshold,"\n")
     return()
  }
  ccorclass1 <- x[[2]][difCorrelated,difCorrelated]
  ccorclass2 <- x[[3]][difCorrelated,difCorrelated]
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
    image(seq(0,nrow(ccorclass1)),seq(0,nrow(ccorclass1)),genotypes,ylab="",xlab="",yaxt="n",xaxt="n",col=c(rgb(255,0,0,50,maxColorValue=255),rgb(255,255,0,50,maxColorValue=255)),main=paste("QCL at marker",attr(x,"marker")),...)
    image(seq(0,nrow(ccorclass1)),seq(0,nrow(ccorclass1)),lower+upper,breaks=c(-1,-0.75,-0.5,0.5,0.75,1),col=c("blue","lightblue",rgb(255,255,255,155,maxColorValue=255),"lightgreen","green"),add=TRUE,...)
    axis(2,at=seq(0.5,nrow(ccorclass1)),labels=colnames(ccorclass1),las=1)
    axis(1,at=seq(0.5,nrow(ccorclass1)),labels=colnames(ccorclass1),las=3)
    if(nrow(genotypes) < 50){
      grid(nrow(ccorclass1),nrow(ccorclass1),lwd=1,lty=1,col="black")
    }else{
      cat("INFO: ",attr(x,"marker"),"- 50+ traits, plot grid disabled, Add with: grid(",nrow(genotypes),",",nrow(genotypes),",lwd=1,lty=1,col=\"black\")\n")
    }
    abline(v=0)
    abline(h=nrow(ccorclass1))
    ordering
  }else{
    warning("No phenotype shows differential correlation with more than: ",QCL.significant," other phenotypes at QCL.threshold: ",QCL.threshold,"\n")
  }
}

#Heatmap the output of a markerQCL object (almost the same as above)
plot.markerQCL2 <- function(x, peekheight=1, ...){
  selection <- which(x[[4]] >= peekheight)
  if(length(selection) > 0){
    g1 <- x[[2]][selection,selection]
    g2 <- x[[3]][selection,selection]
    ordering <- x[[1]][selection,selection]
    clustering <- hclust(dist(ordering))

    upper <- upper.tri(g1)*sign(g1)*(g1[clustering$order,clustering$order]^2)
    lower <- lower.tri(g2)*sign(g2)*(g2[clustering$order,clustering$order]^2)
    colorz <- c("red","yellow","white","lightblue","blue")
    heatmap(t(upper+lower),col=colorz,breaks=c(-1,-0.5,-0.3,0.3,0.5,1),Colv=NA,Rowv=NA,scale="none",main=paste("Correlation at: ",attr(x,"marker")))
    clustering
  }else{
    plot(1:10)
    return(NULL)
  }
}

#Heatmap the output of a QCLscan object
image.QCLscan <- function(x, QCL.thresholds=0, cluster=FALSE, ...){
  if(cluster){
    heatmap(t(x[,which(apply(x,2,function(x){max(x) > QCL.thresholds}))]),Colv=NA, scale="none", col=c("white",gray.colors(10)[10:1],"black"))
  }else{
    heatmap(t(x[,which(apply(x,2,function(x){max(x) > QCL.thresholds}))]),Colv=NA, Rowv=NA, scale="none", col=c("white",gray.colors(10)[10:1],"black"))
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

bplot.QCLscan <- function(x, cross,addQTL=TRUE){
  phenoname <- attr(x,"phenoname")
  offsets <- getMarkerOffsets(cross)
  plot(offsets,as.numeric(c(max(x),rep(min(x),sum(nmar(cross))-1))),type="n",main=phenoname,xlab="",ylab="")
  
  olnm <- 1
  snm <- 0
  if(addQTL){
    ids <- which(phenames(cross) == phenoname)
    cat(phenoname,"=",ids,"\n")
    cross <- calc.genoprob(cross)
    qtlscan <- scanone(cross, pheno.col=ids)
    dtmax <- max(x)
    qtlmax <- max(qtlscan[,3],10)
    qtlscan[,3] <- qtlscan[,3]*(dtmax/qtlmax)
    cat(qtlmax," ",dtmax,"\n")
    axis(4,at=seq(0,1.7*dtmax,1),round((qtlmax/dtmax) * seq(0,1.7*dtmax,1),1))
  }
  for(nm in nmar(cross)){
    snm = nm + snm
    points(offsets[olnm:snm],as.numeric(qtlscan[olnm:snm,3]),lwd=3,type='o',col="red")
    points(offsets[olnm:snm],as.numeric(x[olnm:snm,phenoname]),lwd=2,type='o',col="blue")
    olnm=snm+1
  }
  boxplot(t(x), at=offsets, main=phenoname, add=TRUE, xlab=FALSE, names=rep("",length(offsets)), ylab="QCL",pch=20)
  legend("topright",c("Boxplot QCL","QCL","QTL (scanone)"),lty=c(1,1,1),lwd=c(1,3,2),col=c("black","blue","red"))
}


#plot.QCLscan Plot a QCL profile of a single phenotype
#Optionally add a scanone object to overlay the QTL profile
plot.QCLscan <- function(x, cross, QCL.thresholds=NULL, addQTL=FALSE, ...){
  todrop <- markernames(cross)[unlist(1:sum(nmar(cross)))[-attr(x,"markers")]]
  cross <- drop.markers(cross,todrop)
  phenoname <- attr(x,"phenoname")
  QCLscanProfile <- lodscorestoscanone(cross,apply(x,1,var),traitnames = "QCL")
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
plot.markerQCL.detailed <- function(x, pheno.col=1, qcl.threshold=0, ...){
  sorted <- sort(abs(x[[1]][,pheno.col]),decreasing=TRUE)
  sorted <- sorted[which(sorted > qcl.threshold)]
  ordering <- names(sorted)
  plot(c(1,length(ordering)),c(0,1),xlab="Other phenotypes",ylab="Correlation",main=paste("Correlation probe ",pheno.col," at marker",attr(x,"marker"),sep=" "),type='n')
  points(x[[2]][ordering,pheno.col]^2,pch=20,col='red',type='s')
  points(x[[3]][ordering,pheno.col]^2,pch=20,col='green',type='s')
  points(abs(x[[1]][ordering,pheno.col]),pch=20,col='blue',type='l',lwd=4)
  legend("topright",c("Genotype 1","Genotype 2","Absolute difference"),pch=20,col=c("red","green","blue"))
  invisible(ordering)
}


#colorrange <- c(rgb(seq(0,1,0.01),seq(0,1,0.01),1),"white",rgb(seq(1,0,-0.01),1,seq(1,0,-0.01)))
#for(x in 1:3){
#  sign <- unique(unlist(lapply(apply(QCL[[x]],2,function(y){which(y>0.45)}),names)))
#}
#  plot(apply(QCL[[x]],2,var,na.rm=T),type='l')
#  heatmap(QCL[[x]][sign,],Colv=NA, col=colorrange,scale="none")
#  image(1:ncol(QCL[[x]]),1:nrow(QCL[[x]]),t(QCL[[x]]),col=colorrange)


# Counts the number of ocurences above threshold in a vector
countVThreshold <- function(vector,threshold = 0.5){
  length(which(abs(vector) >= threshold))
}

# Counts the number of occurences above threshold in the QCL
count.QCL <- function(QCL,threshold = 0.5){
  apply(QCL[[1]],1,countVThreshold,threshold)
}

summary.QCL <- function(object, ...){
  for(s in c(seq(5,50,5),100,200,500,1000)){
    not_significant <- apply(object,2,function(x){sum(x>s)==0})
    significant <- apply(object,2,function(x){sum(x>s)!=0})
    cat(s," ",sum(not_significant)," ",sum(significant),"\n")
  }
}

#Test (and time) the differentialCorrelation routine
test.QCL <- function(){
  require(qtl)
  data(multitrait)
  s <- proc.time()
  multitrait = fill.geno(multitrait)
  res <- QCL.scan(multitrait)
  plot(res)
  e <- proc.time()
  
}

#Create a correlation matrix (so we can easily apply
correlationmatrix <- function(x,expressions,method="pearson"){
	cor(expressions[x,], method = method, use="pairwise.complete.obs")
}
