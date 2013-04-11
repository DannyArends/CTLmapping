
getSNP <- function(firstbase){
  if(levels(firstbase)[1]=='0') return(levels(firstbase)[-1])
  return(levels(firstbase))
}

createGenotype_N3 <- function(firstbase,secondbase, x){
  genotype   <- NULL
  snps       <- getSNP(as.factor(firstbase[,x]))
  mallele    <- names(which.max(table(unlist(firstbase[,1],secondbase[,1]))))
  for(b in 1:length(firstbase[,x])){
    if(as.character(firstbase[b,x])==as.character(secondbase[b,x])){
      if(as.character(firstbase[b,x])==snps[1]){
        genotype <- c(genotype,1)
      }else{
        genotype <- c(genotype,3)
      }
    }else{
      genotype <- c(genotype,2)
    }
  }
  return(genotype)
}

createGenotype_N2 <- function(firstbase,secondbase, x){
  genotype   <- NULL
  snps       <- getSNP(as.factor(firstbase[,x]))
  mallele    <- names(which.max(table(unlist(firstbase[,x],secondbase[,x]))))
  for(b in 1:length(firstbase[,x])){
    if(as.character(firstbase[b,x])==as.character(secondbase[b,x])){
      if(as.character(firstbase[b,x])==mallele){
        genotype <- c(genotype,1)
      }else{
        genotype <- c(genotype,NA)
      }
    }else{
      genotype <- c(genotype,2)
    }
  }
  cat(mallele,' ',length(which(genotype==1)),' ',length(which(genotype==2)),'\n')
  return(genotype)
}

harmjanToLude <- function(hjnum = 32268){ koppel[which(koppel[,1]==hjnum),2] }
LudeToHJ <- function(ludenum = 270615){ koppel[which(koppel[,2]==ludenum),1] }

mapQTLs <- function(expdata, firstbase, secondbase){
  qtlmatrix <- NULL
  for(snp in 1:ncol(firstbase)){
    cat('SNP', snp,'\n')
    genotype   <- createGenotype_N3(firstbase, secondbase, snp)
    models     <- aov(apply(expdata,1,as.numeric) ~ as.numeric(as.factor(genotype)))
    modelinfo  <- summary(models)
    qtls       <- -log10(as.numeric(unlist(lapply(modelinfo,"[",1,5))))
    qtlmatrix  <- cbind(qtlmatrix,qtls)
  }
  return(qtlmatrix)
}