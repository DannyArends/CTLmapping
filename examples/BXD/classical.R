#
# Analysis of classical phenotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Sep, 2016
# first written Sep, 2016
#

setwd("D:/Github/CTLmapping/examples/BXD/data")
phenotypes <- t(read.table("classical/phenotypes.txt", sep="\t", row.names=1, header=TRUE))

### Genotypes
genotypes <- read.table("BXD.geno",sep="\t", skip = 6, header =TRUE, row.names = 2, na.strings = c("U", "H"), colClasses="character")
map <- genotypes[,1:3]
genotypes <- genotypes[,-c(1:3)]

renames <- c("BXD96","BXD97","BXD92","BXD80","BXD103")
names(renames) <- c("BXD48a", "BXD65a", "BXD65b", "BXD73a", "BXD73b")

### Match individual names
pii <- which(colnames(phenotypes) %in% names(renames))            # Rename some individuals in the phenotype dataset
colnames(phenotypes)[pii] <- renames[colnames(phenotypes)[pii]]   # Rename some individuals in the phenotype dataset

for(x in 1:(nrow(phenotypes)-1)){
  for(y in x:nrow(phenotypes)){
    hasPhe <- names(which(apply(apply(phenotypes[c(x,y),],2,is.na),2,sum) == 0))
    subPhe <- phenotypes[c(x,y), hasPhe]
    subGen <- genotypes[,colnames(subPhe)]
    
    subGen <- subGen[which(apply(apply(subGen, 1, is.na), 2, sum) == 0),]
    map <- map[rownames(subGen),]
    
    # Change encoding of the genotypes to numeric 1 and 2
    subGen[subGen == "B"] <- 1
    subGen[subGen == "D"] <- 2
    res <- CTLscan(t(subGen), t(subPhe), phenocol=1)
    if(any(apply(res[[1]]$ctl,1,max) > 0)) cat("CTL",x," ",y, " ", max(apply(res[[1]]$ctl,1,max)),"\n")
  }
}