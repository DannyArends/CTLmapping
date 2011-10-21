#
#
# QCL.R
#
# copyright (c) 2010 Danny Arends, Bruno Tesson and Ritsert C. Jansen
# last modified feb, 2011
# first written nov, 2010
# 
# R functions to do QCL mapping

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
  res <- QCLscan(multitrait)
  plot(res)
  e <- proc.time()
  
}

#Create a correlation matrix (so we can easily apply
correlationmatrix <- function(x,expressions,method="pearson"){
	cor(expressions[x,], method = method, use="pairwise.complete.obs")
}

#QCL internal routine - Optimized using snow and 2 cores
#Does a Single marker return a list with:
# [[1]] difCorMatrix 
# [[2]] Correlation matrix individuals with genotype 1
# [[3]] Correlation matrix individuals with genotype 2
# [[4]] A vector of counted scores for each phenotype above the difCorThreshold
#Note: Also saves the object to: output/difCor<marker>.Rdata
QCLscan.internal <- function(cross, marker, QCL.threshold=0.25, method="pearson", directory="output", saveRdata=FALSE, cpu_cluster){
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
QCLscan <- function(cross, pheno.col=1, marker.col, QCL.threshold=0.25, significant = 0, method="pearson", directory="output", doplot=FALSE, writefile=FALSE, saveRdata=FALSE, snow=TRUE, verbose=TRUE){
  s <- proc.time()
  if(doplot && !file.exists(directory)) dir.create(directory)
  if(verbose) cat("Analysis of ",ncol(cross$pheno)," traits at ",sum(nmar(cross))," markers\n")
  phenoname <- phenames(bremcross)[pheno.col]
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
    results <- QCLscan.internal(cross, marker, QCL.threshold, method,directory,saveRdata, cpu_cluster)
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
