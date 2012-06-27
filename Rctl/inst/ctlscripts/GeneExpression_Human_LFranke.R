setwd("E:\\GBIC\\LudeCTL")
memory.limit(2000)

harmjanToLude <- function(hjnum = 32268){ koppel[which(koppel[,1]==hjnum),2] }

expdata    <- read.csv("out_40PCA_noG.txt", sep="\t", row.names=1)
transQTLs  <- read.csv("input/eQTLsFDR0.05.txt", header=TRUE, sep="\t",row.names=NULL)
koppel     <- read.csv("input/koppeltabel1.txt", header=TRUE, sep="\t",row.names=NULL)

hjnames    <- as.character(transQTLs[,5])
allnames   <- unique(as.character(unlist(lapply(hjnames,harmjanToLude))))
expids     <- which(rownames(expdata) %in% allnames)

expdata <- expdata[expids, ]
namesexp <- gsub(".","-",colnames(expdata),fixed=T)
namesexp <- gsub("HT12_","-",namesexp)
colnames(expdata) <- namesexp

getDataMatrix <- function(chr, loc, margin = 100){
  info <- read.csv(paste("genotypes/chr",chr,".map",sep=""), header=FALSE, sep="\t", row.names=2)
  ids <- which(loc-margin < info[,3] & info[,3] < loc+margin)
  cat("Going to select:", length(ids)," SNP markers\n")
  genotype_info <- info[ids,]
  ids <- ((min(ids)*2)-1):(max(ids)*2)
  genotypes <- read.csv(paste("genotypes/chr",chr,".ped",sep=""), header=FALSE, sep=" ", row.names=2)[,-c(1:5)]
  list(genotype_info,genotypes[,ids])
}

chr2data   <- getDataMatrix(2,   60573474, 250000)
chr3data   <- getDataMatrix(3,   56840816, 250000)
chr6_1data   <- getDataMatrix(6,  139880112, 250000)
chr6_2data   <- getDataMatrix(6,  135468837, 250000)
chr7data   <- getDataMatrix(7,   50395922, 250000)
chr11data   <- getDataMatrix(11,    192856, 250000)
chr12_1data   <- getDataMatrix(12,  54756892, 250000)
chr12_2data   <- getDataMatrix(12, 110368991, 250000)

red_genotypes <- cbind(chr2data[[2]],chr3data[[2]],chr6_1data[[2]],chr6_2data[[2]],chr7data[[2]],chr11data[[2]],chr12_1data[[2]],chr12_2data[[2]])
genonames <- rbind(chr2data[[1]][,c(1,3)],chr3data[[1]][,c(1,3)],chr6_1data[[1]][,c(1,3)],chr6_2data[[1]][,c(1,3)],chr7data[[1]][,c(1,3)],chr11data[[1]][,c(1,3)],chr12_1data[[1]][,c(1,3)],chr12_2data[[1]][,c(1,3)])

firstbase     <- red_genotypes[, seq(1,ncol(red_genotypes),2)]
secondbase    <- red_genotypes[, seq(2,ncol(red_genotypes),2)]

ordering <- match(rownames(chr12_2data[[2]]), colnames(expdata))
expdata <- expdata[,ordering]

#Check if the colnames of the expression data match the colnames of the rawgenotypes
colnames(expdata)[1:50]
rownames(chr12_2data[[2]])[1:50]

getSNP <- function(firstbase){
  if(levels(firstbase)[1]=='0') return(levels(firstbase)[-1])
  return(levels(firstbase))
}

createGenotype <- function(firstbase,secondbase, x){
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

genotypes <- NULL
qtlmatrix <- NULL
for(x in 1:ncol(firstbase)){
  genotype   <- createGenotype(firstbase,secondbase,x)
 # models     <- aov(apply(expdata,1,as.numeric) ~ as.numeric(as.factor(genotype)))
  #modelinfo  <- summary(models)
  #qtls       <- -log10(as.numeric(unlist(lapply(modelinfo,"[",1,5))))
  #qtlmatrix  <- cbind(qtlmatrix,qtls)
  genotypes  <- cbind(genotypes,genotype)

 # plot(c(0,nrow(expdata)),c(0,50), t='n', xlab="Probe", ylab="LOD")
 # points(qtls, t='h', xlab="Probe", ylab="LOD")
}

image(t(qtlmatrix),breaks=c(0,5,10,100),col=c("white","gray","black"),xlab="Loc (bp)",ylab="Probe")
box()

sums <- apply(genotypes,2,function(x){c(sum(x==3,na.rm=T),sum(x==1,na.rm=T))})
bad_markers <- unique(c(which(sums[1,] < 20),which(sums[2,] < 20)))

write.table(file="n3genotypes.txt",cbind(genonames[-bad_markers,],t(genotypes[,-bad_markers])),sep="\t",na="",col.names=FALSE)
write.table(file="n3phenotypes.txt",cbind(NA,NA,expdata),sep="\t",na="",col.names=FALSE,quote=FALSE)

