#
# Analysis of classical phenotype data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Sep, 2016
# first written Sep, 2016
#
library(ctl)

setwd("D:/Github/CTLmapping/examples/BXD/data")
phenotypes <- read.csv("classical/export-16-09-13-06-34.txt", sep="\t", header=TRUE, skip=10,na.strings=c("", "None"))

idx <- grep("_Value", colnames(phenotypes))
colnames(phenotypes)[idx] <- gsub("_Value", "", colnames(phenotypes)[idx])

pubmedid <- phenotypes[,10]
phenotypes <- phenotypes[,idx[-c(1:4)]]

renames <- c("BXD96","BXD97","BXD92","BXD80","BXD103")
names(renames) <- c("BXD48a", "BXD65a", "BXD65b", "BXD73a", "BXD73b")

### Match individual names
pii <- which(colnames(phenotypes) %in% names(renames))            # Rename some individuals in the phenotype dataset
colnames(phenotypes)[pii] <- renames[colnames(phenotypes)[pii]]   # Rename some individuals in the phenotype dataset


### Genotypes
genotypes <- read.table("BXD.geno",sep="\t", skip = 6, header =TRUE, row.names = 2, na.strings = c("U", "H"), colClasses="character")
map <- genotypes[,1:3]
genotypes <- genotypes[,-c(1:3)]
interesting <- NULL

for(x in 1:(nrow(phenotypes)-1)){
  for(y in 56:56){
    hasPhe <- names(which(apply(apply(phenotypes[c(x,y),],2,is.na),2,sum) == 0))
    subPhe <- phenotypes[c(x,y), hasPhe]
    subGen <- genotypes[,colnames(subPhe)]

    subGen <- subGen[which(apply(apply(subGen, 1, is.na), 2, sum) == 0),]             # No missing genotypes
    subGen <- subGen[which(lapply(apply(subGen,1,table), min) > 5),]                  # At least 5 indidivudals in each group
    map <- map[rownames(subGen),]
    err <- 0
    # Change encoding of the genotypes to numeric 1 and 2
    subGen[subGen == "B"] <- 1
    subGen[subGen == "D"] <- 2
    tryCatch(res <- CTLscan(t(subGen), t(subPhe), phenocol=1), error = function(e){err <- 1})
    if(err == 0 && any(apply(res[[1]]$ctl,1,max) > 0)){
      interesting <- cbind(interesting,c(x,y))
      cat("CTL", x," ", y, " ", max(apply(res[[1]]$ctl,1,max)),"\n")
    }
    cat(y,"\n")
  }
}

