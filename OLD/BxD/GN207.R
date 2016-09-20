#
# Analysis of GN207 BxD data
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

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

### Phenotypes
BXDdata <- read.csv(gzfile("GN207/GN207_MeanDataAnnotated_rev081815.txt.gz"), skip = 33, header = TRUE, sep="\t", colClasses="character")
phenotypes <- BXDdata[which(BXDdata[,"Gene.Symbol"] %in% highImpact[,1]),]
annotation <- phenotypes[,c("Gene.Symbol", "Description")]

### Match individual names
pii <- which(colnames(phenotypes) %in% names(renames))            # Rename some individuals in the phenotype dataset
colnames(phenotypes)[pii] <- renames[colnames(phenotypes)[pii]]   # Rename some individuals in the phenotype dataset

# Subset phenotypes and genotypes
phenotypes  <- phenotypes[, which(colnames(phenotypes) %in% colnames(genotypes))]
genotypes   <- genotypes[, colnames(phenotypes)]

# Complete genotype data
genotypes <- genotypes[which(apply(apply(genotypes,1,is.na),2,sum) == 0),]
map <- map[rownames(genotypes),]

# Change encoding of the genotypes to numeric 1 and 2
genotypes[genotypes == "B"] <- 1
genotypes[genotypes == "D"] <- 2

phenotypes[1:5, 1:15]
annotation[1:5, ]
genotypes[1:5, 1:15]

# Significance levels
qtl_cutoff <- -log10(0.05 / nrow(genotypes))
ctl_cutoff <- -log10(0.05)

# Map QTLs
if(!file.exists("GN207/QTLs.txt")) {
  QTLs <- matrix(NA, nrow(phenotypes), nrow(genotypes))
  colnames(QTLs) <- rownames(genotypes)

  for(x in rownames(genotypes)){
    aovModel <- aov(apply(phenotypes, 1, as.numeric) ~ as.character(genotypes[x,]))
    QTLs[,x] <- unlist(lapply(summary(aovModel),function(x){ return(x[[5]][1]) }))
  }

  write.table(QTLs, file = "GN207/QTLs.txt", sep="\t", quote=FALSE)
}else{
  QTLs <- read.table("GN207/QTLs.txt")
}

# Map CTL, here we can only use a subset at a time ( I do 100, then close R and continue)
if(!file.exists("GN207/CTLs_p.txt")) {
  CTLs <- matrix(NA, nrow(phenotypes), nrow(genotypes))
  colnames(CTLs) <- rownames(genotypes)

  cat("", file="GN207/CTLs_p.txt")
  cat("", file="GN207/CTLs_int.txt")
  for(x in 1:nrow(phenotypes)){
    res <- CTLscan(t(genotypes), t(phenotypes), phenocol=x)                                             # Scan for CTLs
    scores <- apply(res[[1]]$ctl, 1, max)                                                               # Max CTL scores per marker
    significant <- which(apply(res[[1]]$ctl, 2, max) > ctl_cutoff)                                      # Other phenotypes causing the CTL
    for(y in significant){
      maxM <- which.max(res[[1]]$ctl[,y])
      crs <- getCorrelations(t(genotypes), t(phenotypes), x, maxM, verbose=FALSE)
      gphes <- names(which(res[[1]]$ctl[,y] > ctl_cutoff))                                              # Markers causing the CTL
      locs  <- unique(cbind(map[gphes,1], round(as.numeric(map[gphes,3])/5) * 5))                       # Approx location
      cors <- round(as.numeric(crs$correlations[y,]),2)
      ssize <- as.numeric(crs$samplesize)
      interacts <- cbind(annotation[x,1], annotation[y,1], locs, 
                         round(max(res[[1]]$ctl[gphes,y]),2), 
                         round(mean(res[[1]]$ctl[gphes,y]),2),
                         length(gphes),
                         names(maxM),cors[1],cors[2], ssize[1], ssize[2], "\n")
      cat(t(interacts), file="GN207/CTLs_int.txt", append = TRUE)
    }
    cat(paste0(paste0(c(x, scores),collapse="\t"), "\n"), file="GN207/CTLs_p.txt", append = TRUE)
    CTLs[x,] <- scores
    if(x %% 4 == 0) cat(x, "\n")
  }
  write.table(CTLs, file = "GN207/CTLs.txt", sep="\t", quote=FALSE)
}else{
  CTLs <- read.table("GN207/CTLs_p.txt", row.names=1)
}

haveQTL <- which(apply(-log10(QTLs), 1, max) > qtl_cutoff)           # 5 is 'too low' for QTL
haveCTL <- which(apply(CTLs, 1, max) > ctl_cutoff)                   # 5 is 'too high/stringent' for CTL

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

table(unlist(lapply(annotation[,1], whichGroup, highImpact)))
table(unlist(lapply(annotation[haveQTL,1], whichGroup, highImpact)))
table(unlist(lapply(annotation[haveCTL,1], whichGroup, highImpact)))

# From http://stackoverflow.com/questions/2261079
trim.trailing <- function (x) sub("\\s+$", "", x) # returns string w/o trailing whitespace
trim.leading <- function (x)  sub("^\\s+", "", x) # returns string w/o leading whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)   # returns string w/o leading or trailing whitespace

# Create output for excel
# Read in the interaction table and format for cytoscape
interactions <- trim.trailing(readLines("GN207/CTLs_int.txt"))

sourc <- unlist(lapply(strsplit(interactions, " "),"[",1))
targe <- unlist(lapply(strsplit(interactions, " "),"[",2))
chr <- unlist(lapply(strsplit(interactions, " "),"[",3))
pos <- unlist(lapply(strsplit(interactions, " "),"[",4))

uniques <- which(!duplicated(apply(apply(cbind(sourc,targe,chr,pos),1,sort),2,paste0,collapse="-")))
cat(gsub(" ","\t",trim(readLines("GN207/CTLs_int.txt")[uniques])),sep="\n", file="GN207/CTLs_int_all.txt")

# Read in the interaction table and format for cytoscape
interactions <- trim.trailing(readLines("GN207/CTLs_int.txt")[which(substr(readLines("GN207/CTLs_int.txt"),1,1) != " ")])

cat(gsub(" ","\t",trim(readLines("GN207/CTLs_int.txt"))),sep="\n", file="GN207/CTLs_int_all.txt")

sourc <- unlist(lapply(strsplit(interactions, " "),"[",1))
targe <- unlist(lapply(strsplit(interactions, " "),"[",2))
chr <- unlist(lapply(strsplit(interactions, " "),"[",3))
pos <- unlist(lapply(strsplit(interactions, " "),"[",4))

uniques <- which(!duplicated(apply(apply(cbind(sourc,targe,chr,pos),1,sort),2,paste0,collapse="-")))
cat(gsub(" ","\t",trim(readLines("GN207/CTLs_int.txt")[uniques])),sep="\n", file="GN207/CTLs_int_unique.txt")

groupS <- unlist(lapply(sourc, whichGroup, highImpact))
groupT <- unlist(lapply(targe, whichGroup, highImpact))
interactions  <- paste0(interactions," ",paste0(chr,":",pos), " ", groupS, " ", groupT)

edges <- NULL
for(x in 1:length(sourc)){
  if(sourc[x] != targe[x])edges <- c(edges, sourc[x], targe[x])
}
g1 <- graph( edges=edges, directed=F )
g1 <- simplify(g1, remove.multiple = T) 
l <- layout_with_mds(g1)
plot(g1,layout=l)

ctlsign <- which(as.numeric(unlist(lapply(strsplit(interactions, " "),"[",5))) > ctl_cutoff)
cat(gsub(" ","\t", c("Source\tTarget\tChr\tPos\tmaxLOD\tmeanLOD\tnSNPs\ttop_marker\tcor1\tcor2\tss1\tss2\tChrPos\tGroup\tGroup", interactions[ctlsign])), sep="\n", file="GN207/CTLs_cyto.txt")

### Figure 1 for publication in JOSS

png("Fig1.png", width=2048, height=700, res=300, pointsize = 5)

op <- par(mfrow=c(1,3),cex = 0.9, mgp = c(2, 1, 0), lwd=0.9)

ctlpr <- as.numeric(CTLs[which(annotation[,1] == "St7")[1],])
qtlpr1 <- -log10(as.numeric(QTLs[which(annotation[,1] == "St7")[1],]))
qtlpr2 <- -log10(as.numeric(QTLs[which(annotation[,1] == "Il18r1")[1],]))
onmap <- which(map[,1] == map[which.max(ctlpr),1])
chrmap <- map[onmap,]
x <- as.numeric(chrmap[,2])

plot(c(0,max(x)), c(-5, 5), t ='n', main = "a) CTL, No QTL : St7 vs Il18r1", xlab="Chromosome 2 - Position (Mb)", ylab = "<- CTL   -log10(P-value)   QTL ->")
points(x, -ctlpr[onmap],t = 'l', col="steelblue")
points(x, -ctlpr[onmap], pch=10, col="steelblue")
points(x, qtlpr1[onmap], t = 'l', lty=1)
points(x, qtlpr1[onmap], pch=19)
points(x, qtlpr2[onmap], t = 'l', lty=2, col="orange")
points(x, qtlpr2[onmap], pch=1, col="orange")
abline(h = -ctl_cutoff, lty=3)
legend("topleft", c("QTL St7", "QTL Il18r1", "CTL St7 vs Il18r1", "5% FDR"), pch=c(19,1,10,NA), lwd=0.9, col=c(1,"orange","steelblue",1), lty=c(1, 2, 1, 3))


ctlpr <- as.numeric(CTLs[which(annotation[,1] == "St7")[1],])
qtlpr1 <- -log10(as.numeric(QTLs[which(annotation[,1] == "St7")[1],]))
qtlpr2 <- -log10(as.numeric(QTLs[which(annotation[,1] == "Il18r1")[2],]))

onmap <- which(map[,1] == map[which.max(qtlpr1),1])
chrmap <- map[onmap,]
x <- as.numeric(chrmap[,2])

plot(c(0,max(x)), c(-12, 12), t ="n", main = "b) QTL, No CTL : St7 vs Il18r1", xlab="Chromosome 6 - Position (Mb)", ylab = "<- CTL   -log10(P-value)   QTL ->")
points(x, -ctlpr[onmap],t = 'l', col="steelblue")
points(x, -ctlpr[onmap], pch=10, col="steelblue")
points(x, qtlpr1[onmap], t = 'l', lty=1)
points(x, qtlpr1[onmap], pch=19)
points(x, qtlpr2[onmap], t = 'l', lty=2, col="orange")
points(x, qtlpr2[onmap], pch=1, col="orange")
abline(h = qtl_cutoff, lty=3)
abline(h = -ctl_cutoff, lty=3)
legend("topright", c("QTL St7", "QTL Il18r1", "CTL St7 vs Il18r1", "5% FDR"), pch = c(19,1,10,NA), col=c(1,"orange","steelblue",1), lwd=0.9, lty=c(1, 2, 1, 3))


ctlpr <- as.numeric(CTLs[which(annotation[,1] == "Mtvr2")[1],])
qtlpr1 <- -log10(as.numeric(QTLs[which(annotation[,1] == "Mtvr2")[1],]))
qtlpr2 <- -log10(as.numeric(QTLs[which(annotation[,1] == "C1qtnf5")[1],]))

onmap <- which(map[,1] == map[which.max(ctlpr),1])
chrmap <- map[onmap,]
x <- as.numeric(chrmap[,2])

plot(c(0,max(x)), c(-12, 12), t ='n', main = "c) QTL and CTL : Mtvr2 vs C1qtnf5", xlab="Chromosome 19 - Position (Mb)", ylab = "<- CTL   -log10(P-value)   QTL ->")
points(x, -ctlpr[onmap],t = 'l', col="steelblue")
points(x, -ctlpr[onmap], pch=10, col="steelblue")
points(x, qtlpr1[onmap], t = 'l', lty=1)
points(x, qtlpr1[onmap], pch=19)
points(x, qtlpr2[onmap], t = 'l', lty=2, col="orange")
points(x, qtlpr2[onmap], pch=1, col="orange")
abline(h = qtl_cutoff, lty=3)
abline(h = -ctl_cutoff, lty=3)
legend("topright", c("QTL Mtvr2", "QTL C1qtnf5", "CTL Mtvr2 vs C1qtnf5", "5% FDR"), pch = c(19,1,10,NA), col=c(1,"orange","steelblue",1), lwd=0.9, lty=c(1, 2, 1, 3))

dev.off()
