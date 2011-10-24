
QCLscan <- function(genotypes,phenotypes, pheno.col, verbose = F){
  cnt <- 1
  results <- vector("list",length(pheno.col))
  for(x in pheno.col){
    profile <- apply(genotypes,2, 
      function(geno){
        cor1 <- cor(phenotypes[geno==1,x],phenotypes[geno==1,],use="pair")
        cor2 <- cor(phenotypes[geno==2,x],phenotypes[geno==2,],use="pair")
        sign(cor1)*(cor1^2)-sign(cor2)*(cor2^2)
      }
    )
    rownames(profile) <- colnames(phenotypes)
    colnames(profile) <- colnames(genotypes)
    results[[cnt]] <- profile
    attr(results[[cnt]],"name") <- colnames(phenotypes)[x]
    if(verbose){
      cat("Phenotype:",colnames(phenotypes)[x],"\n")
    }
    cnt <- cnt +1
  }
  results
}

QCLscanToSIF <- function(QCLscan, cutoff=0.45){
  cat("",file="network.sif")
  cnt <- 0
  for(QCL in QCLscan){
    for(x in 1:ncol(QCL)){
      for(id in which(QCL[,x] > cutoff)){
        cat(attr(QCL,"name"),"QCL",rownames(QCL)[id],QCL[id,x],colnames(QCL)[x],"\n",file="network.sif",append=TRUE)
        cnt <- cnt+1
      }
    }
    cat("Next phenotype\n")
  }
  cat("Wrote",cnt,"Edges\n")
}

cor1 <- c(-1,-1,-1, 0, 0, 0, 1, 1, 1)
cor2 <- c(-1, 0, 1,-1, 0, 1,-1, 0, 1)

phenotypes <- apply(pull.pheno(cross),2,as.numeric)
genotypes <- pull.geno(cross)

QCL <- QCLscan(genotypes,phenotypes,1:20,verbose=T)

colorrange <- c(rgb(seq(0,1,0.01),seq(0,1,0.01),1),"white",rgb(seq(1,0,-0.01),1,seq(1,0,-0.01)))


for(x in 1:3){
  sign <- unique(unlist(lapply(apply(QCL[[x]],2,function(y){which(y>0.45)}),names)))
#  plot(apply(QCL[[x]],2,var,na.rm=T),type='l')
  heatmap(QCL[[x]][sign,],Colv=NA, col=colorrange,scale="none")
  image(1:ncol(QCL[[x]]),1:nrow(QCL[[x]]),t(QCL[[x]]),col=colorrange)
}

setwd("e:/")
QCLscanToSIF(QCL)