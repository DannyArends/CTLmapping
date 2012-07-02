setwd("E:\\GBIC\\Lude CTL")
memory.limit(2000)

harmjanToLude <- function(hjnum = 32268){ koppel[which(koppel[,1]==hjnum),2] }

expdata    <- read.csv("out_40PCA_noG_eQTL.txt", sep="\t", row.names=1,check.names = FALSE)
transQTLs  <- read.csv("input/eQTLsFDR0.05.txt", header=TRUE, sep="\t",row.names=NULL)
koppel     <- read.csv("input/koppeltabel1.txt", header=TRUE, sep="\t",row.names=NULL)

hjnames    <- unique(as.character(transQTLs[,5]))
allnames   <- unique(as.character(unlist(lapply(hjnames,harmjanToLude))))

getDataMatrix <- function(chr, loc, margin = 100){
  info <- read.csv(paste("genotypes/chr",chr,".map",sep=""), header=FALSE, sep="\t", row.names=2)
  ids <- which(loc-margin < info[,3] & info[,3] < loc+margin)
  cat("Going to select:", length(ids)," SNP markers\n")
  genotype_info <- info[ids,]
  ids <- ((min(ids)*2)-1):(max(ids)*2)
  genotypes <- read.csv(paste("genotypes/chr",chr,".ped",sep=""), header=FALSE, sep=" ", row.names=2)[,-c(1:5)]
  list(genotype_info,genotypes[,ids])
}

getSNP <- function(firstbase){
  if(levels(firstbase)[1]=='0') return(levels(firstbase)[-1])
  return(levels(firstbase))
}

createGenotype_N3 <- function(firstbase,secondbase, x){
  genotype   <- NULL
  snps       <- getSNP(firstbase[,x])
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

table(apply(transQTLs[,c(3,4)],1,function(x){paste(x,collapse=" ")}))

chr3data   <- getDataMatrix(12, 54768447, 2000)

#Match to the genotypes
expdata <-  expdata[, match(rownames(chr3data[[2]]), colnames(expdata))]
#Match to the transQTLs
expdata <-  expdata[match(allnames,rownames(expdata)),]

qtlmatrix <- NULL
firstbase     <- chr3data[[2]][, seq(1,ncol(chr3data[[2]]),2)]
secondbase    <- chr3data[[2]][, seq(2,ncol(chr3data[[2]]),2)]
if(is.null(ncol(firstbase))){
  firstbase <- t(t(firstbase))
  secondbase <- t(t(secondbase))
}
for(snp in 1:ncol(firstbase)){
  cat('SNP', snp,'\n')
  genotype   <- createGenotype_N3(firstbase,secondbase,snp)
  models     <- aov(apply(expdata,1,as.numeric) ~ as.numeric(as.factor(genotype)))
  modelinfo  <- summary(models)
  qtls       <- -log10(as.numeric(unlist(lapply(modelinfo,"[",1,5))))
  qtlmatrix  <- cbind(qtlmatrix,qtls)
}

image(qtlmatrix,breaks=c(0,20,1000),col=c('white','black'))
box()
max(qtlmatrix)
