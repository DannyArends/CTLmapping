#
# Analysis of GN206 BxD data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

### Genotypes
setwd("D:/Github/CTLmapping/examples/BXD/data")

genotypes <- read.table("BXD.geno",sep="\t", skip = 6, header =TRUE, row.names = 2, na.strings = c("U", "H"), colClasses="character")
map <- genotypes[,1:3]
genotypes <- genotypes[,-c(1:3)]

renames <- c("BXD96","BXD97","BXD92","BXD80","BXD103")
names(renames) <- c("BXD48a", "BXD65a", "BXD65b", "BXD73a", "BXD73b")

BXDdata <- read.csv(gzfile("GN206/GN206_MeanDataAnnotated_rev081815.txt.gz"), skip = 33, header = TRUE, sep="\t", colClasses="character")
BXDdata[1:10,]

highImpact <- readLines("genes.txt")


phenotypes <- BXDdata[which(BXDdata[,"Gene.Symbol"] %in% highImpact),]
annotation <- phenotypes[,c("Gene.Symbol", "Description")]
pii <- which(colnames(phenotypes) %in% names(renames))            # Rename some individuals in the phenotype dataset
colnames(phenotypes)[pii] <- renames[colnames(phenotypes)[pii]]   # Rename some individuals in the phenotype dataset

# Subset phenotypes and genotypes
phenotypes  <- phenotypes[, which(colnames(phenotypes) %in% colnames(genotypes))]
genotypes   <- genotypes[, colnames(phenotypes)]

phenotypes[1:5, 1:15]
annotation[1:5, ]
genotypes[1:5, 1:15]

# Map QTLs
if(!file.exists("GN206/QTLs.txt")) {
  QTLs <- matrix(NA, nrow(phenotypes), nrow(genotypes))
  colnames(QTLs) <- rownames(genotypes)

  for(x in rownames(genotypes)){
    aovModel <- aov(apply(phenotypes, 1, as.numeric) ~ as.character(genotypes[x,]))
    QTLs[,x] <- unlist(lapply(summary(aovModel),function(x){ return(x[[5]][1]) }))
  }

  write.table(QTLs, file = "GN206/QTLs.txt", sep="\t", quote=FALSE)
}else{
  QTLs <- read.table("GN206/QTLs.txt")
}

# Write out the data as a cross object for R/qtl
phenotypes <- rbind(phenotypes, sex=rep("M", ncol(phenotypes)))
write.table(cbind(NA,NA, phenotypes), "GN206/cross.csvr", sep=",", quote=FALSE, col.names=FALSE, row.names=TRUE, na="")
write.table(cbind(map[,1],map[,2], genotypes), "GN206/cross.csvr", sep=",", append=TRUE, quote=FALSE, col.names=FALSE, row.names=TRUE, na="")
library(qtl)
cross <- read.cross("csvr", file="GN206/cross.csvr")

# Write out the data as a cross object for R/qtl
library(ctl)
BXDdata <- NULL
cross <- fill.geno(cross)
gc(TRUE)
# Map CTL, here we can only use a subset at a time ( I do 100, then close R and continue)
if(!file.exists("GN206/CTLs_p.txt")) {
  CTLs <- matrix(NA, nrow(phenotypes), nrow(genotypes))
  colnames(CTLs) <- rownames(genotypes)

  cat("", file="GN206/CTLs_p.txt")
  cat("", file="GN206/CTLs_int.txt")
  for(x in 1:100){
    res   <- CTLscan.cross(cross, phenocol = x)                                               # Scan for CTLs
    scores <- apply(res[[1]]$ctl, 1, max)                                                     # Max CTL scores per marker
    significant <- which(apply(res[[1]]$ctl, 2, max) > 2)                                     # Other phenotypes causing the CTL
    for(y in significant){
      gphes <- names(which(res[[1]]$ctl[,y] > 2))                                             # Markers causing the CTL
      locs  <- unique(cbind(map[gphes,1], round(as.numeric(map[gphes,3])/5) * 5))             # Approx location
      interacts <- cbind(annotation[x,1],annotation[y,1],locs, mean(res[[1]]$ctl[gphes,y]), "\n")
      cat(t(interacts), file="GN206/CTLs_int.txt", append = TRUE)
    }
    cat(paste0(paste0(c(x, scores),collapse="\t"), "\n"), file="GN206/CTLs_p.txt", append = TRUE)
    CTLs[x,] <- scores
    if(x %% 4 == 0) cat(x, "\n")
    gc(TRUE)
  }
  write.table(CTLs, file = "GN206/CTLs.txt", sep="\t", quote=FALSE)
}else{
  CTLs <- read.table("GN206/CTLs_p.txt", row.names=1)
}

haveQTL <- which(apply(-log10(QTLs), 1, max) > 5)           # 5 is 'too low' for QTL
haveCTL <- which(apply(CTLs, 1, max) > 5)                   # 5 is 'too high/stringent' for CTL

dev.off()   # Close any open device and start plotting the CTL profile per probe
for(x in unique(haveQTL, haveCTL)) {
  png(paste0("img/plot",x,".png"))
  plot(c(0,ncol(QTLs)), c(-5, 10), t = 'n', main=annotation[x,1])
  points(-log10(as.numeric(QTLs[x,])), t = 'h', col=as.numeric(as.factor(map[,1])))
  abline(h=4, lty=2)
  points(-as.numeric(CTLs[x,]), t = 'h', col=as.numeric(as.factor(map[,1])))
  abline(h=-2, lty=2)
  dev.off()
}

for(x in haveQTL) {
  cat(annotation[x,1], whichGroup(annotation[x,1]), max(-log10(QTLs[x,])), "\n")
}

for(x in haveCTL) {
  cat(annotation[x,1], whichGroup(annotation[x,1]), max(CTLs[x,]), "\n")
}

# From http://stackoverflow.com/questions/2261079
trim.trailing <- function (x) sub("\\s+$", "", x) # returns string w/o trailing whitespace
trim.leading <- function (x)  sub("^\\s+", "", x) # returns string w/o leading whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)   # returns string w/o leading or trailing whitespace

# Create output for excel
# Read in the interaction table and format for cytoscape
interactions <- trim.trailing(readLines("GN206/CTLs_int.txt"))

sourc <- unlist(lapply(strsplit(interactions, " "),"[",1))
targe <- unlist(lapply(strsplit(interactions, " "),"[",2))
chr <- unlist(lapply(strsplit(interactions, " "),"[",3))
pos <- unlist(lapply(strsplit(interactions, " "),"[",4))

uniques <- which(!duplicated(apply(apply(cbind(sourc,targe,chr,pos),1,sort),2,paste0,collapse="-")))
cat(gsub(" ","\t",trim(readLines("GN206/CTLs_int.txt")[uniques])),sep="\n", file="GN206/CTLs_int_all.txt")

# Read in the interaction table and format for cytoscape
interactions <- trim.trailing(readLines("GN206/CTLs_int.txt")[which(substr(readLines("GN206/CTLs_int.txt"),1,1) != " ")])

sourc <- unlist(lapply(strsplit(interactions, " "),"[",1))
targe <- unlist(lapply(strsplit(interactions, " "),"[",2))
chr <- unlist(lapply(strsplit(interactions, " "),"[",3))
pos <- unlist(lapply(strsplit(interactions, " "),"[",4))

uniques <- which(!duplicated(apply(apply(cbind(sourc,targe,chr,pos),1,sort),2,paste0,collapse="-")))
cat(gsub(" ","\t",trim(readLines("GN206/CTLs_int.txt")[uniques])),sep="\n", file="GN206/CTLs_int_all.txt")

groupS <- unlist(lapply(sourc, whichGroup))
groupT <- unlist(lapply(targe, whichGroup))
interactions  <- paste0(interactions," ",paste0(chr,":",pos), " ", groupS, " ", groupT)




a5 <- which(as.numeric(unlist(lapply(strsplit(interactions, " "),"[",5))) > 5)

cat(gsub(" ","\t", c("Source\tTarget\tChr\tPos\tStrength\tChrPos\tGroup\tGroup", interactions[a5])), sep="\n", file="GN206/CTLs_cyto.txt")
