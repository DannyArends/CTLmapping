
### Setup
setwd("D:/Github/CTLmapping/examples/BXD/data")         # Where is the data ?
source("../whichGroup.R")
library(ctl)

### Selected genes
highImpact <- read.table("genes.txt", sep="\t")

### Genotypes
genotypes <- read.table("BXD.geno",sep="\t", skip = 6, header =TRUE, row.names = 2, na.strings = c("U", "H"), colClasses="character")
map <- genotypes[,1:3]
genotypes <- genotypes[,-c(1:3)]

renames <- c("BXD96","BXD97","BXD92","BXD80","BXD103")
names(renames) <- c("BXD48a", "BXD65a", "BXD65b", "BXD73a", "BXD73b")

#### LOAD Classical phenotypes

setwd("D:/Github/CTLmapping/examples/BXD/data")
phenotypes <- read.csv("classical/export-16-09-13-06-34.txt", sep="\t", header=TRUE, skip=10,na.strings=c("", "None"))

idx <- grep("_Value", colnames(phenotypes))
colnames(phenotypes)[idx] <- gsub("_Value", "", colnames(phenotypes)[idx])

pubmedid <- phenotypes[,c(1:10)]
phenotypes <- phenotypes[,idx[-c(1:4)]]

renames <- c("BXD96","BXD97","BXD92","BXD80","BXD103")
names(renames) <- c("BXD48a", "BXD65a", "BXD65b", "BXD73a", "BXD73b")

### Match individual names
pii <- which(colnames(phenotypes) %in% names(renames))            # Rename some individuals in the phenotype dataset
colnames(phenotypes)[pii] <- renames[colnames(phenotypes)[pii]]   # Rename some individuals in the phenotype dataset

#####

#Find good phenotypes to use with this dataset

gooddata <- which(apply(phenotypes[,colnames(genotypes)],1,function(x){sum(!is.na(x))}) > 70)

phenotypes <- phenotypes[gooddata,]
pubmedid <- pubmedid[gooddata,]

vardata <- apply(phenotypes, 1, function(x){var(x,na.rm=TRUE) > 5})

phenotypes <- phenotypes[vardata,]
pubmedid <- pubmedid[vardata,]

hasPhe <- names(which(apply(apply(phenotypes,2,is.na),2,sum) == 0))
subPhe <- phenotypes[, hasPhe]
subGen <- genotypes[,colnames(subPhe)]
subGen[subGen == "B"] <- 1
subGen[subGen == "D"] <- 2

res <- CTLscan(t(subGen), t(subPhe), verbose=TRUE)


for(x in 1:(nrow(phenotypes)-1)){
  for(y in (x+1):nrow(phenotypes)){
    err <- 0
    # Change encoding of the genotypes to numeric 1 and 2
    tryCatch(, error = function(e){err <- 1})
    if(err == 0 && any(apply(res[[1]]$ctl,1,max) > 0)){
      cat("CTL", x," ", y, " ", max(apply(res[[1]]$ctl,1,max)),"\n")
    }
    cat(y,"\n")
  }
}

