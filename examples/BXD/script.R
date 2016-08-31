#
# Analysis of BxD data
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

GN379 <- read.csv(gzfile("GN111_MeanDataAnnotated_rev081815.txt.gz"), skip = 33, header = TRUE, sep="\t", colClasses="character")
GN379[1:10,]

GN379 <- GN379[- which(grepl("AFFX", GN379[,"Gene.Symbol"])), ]
GN379 <- GN379[- which(grepl("Rik", GN379[,"Gene.Symbol"])), ]
GN379 <- GN379[- which(grepl("Gm", GN379[,"Gene.Symbol"])), ]

# CD <- CD4/CD8 antigen // 126 genes
CD <- GN379[c(which(grepl("Cd", GN379[,"Gene.Symbol"]) & 
              grepl("antigen", GN379[,"Description"]) & 
             !grepl("-like", GN379[,"Description"]) & 
             !grepl("cyclin", GN379[,"Description"]) & 
             !grepl("RIKEN", GN379[,"Description"]) & 
             !grepl("cadherin", GN379[,"Description"]) & 
             !grepl("homolog", GN379[,"Description"]) & 
             !grepl("division", GN379[,"Description"]) & 
             !grepl("family", GN379[,"Description"]))), c("Gene.Symbol")]
# IL <- Interleukines
IL <- GN379[c(which(grepl("interleukin ", GN379[,"Description"]) & 
             !grepl("-like", GN379[,"Description"]) & 
             !grepl("-containing", GN379[,"Description"]) & 
             !grepl("binding protein", GN379[,"Description"]) & 
             !grepl("accessory", GN379[,"Description"]) & 
             !grepl("antagonist", GN379[,"Description"]) & 
             !grepl("transducer", GN379[,"Description"]) & 
             !grepl(", alpha", GN379[,"Description"]) & 
             !grepl(", beta", GN379[,"Description"]) & 
             !grepl("enhancer", GN379[,"Description"]) & 
             !grepl("family", GN379[,"Description"]))), c("Gene.Symbol")]
# HLA <- Major histocompatibility complex (MHC / H2)
H2 <- GN379[c(which(grepl("H2-", GN379[,"Gene.Symbol"]))), c("Gene.Symbol")]
# ACE <-  Angiotensin-converting enzyme
ACE <- GN379[c(which(GN379[,"Gene.Symbol"] == "Ace"),
               which(GN379[,"Gene.Symbol"] == "Ace2"),
               which(GN379[,"Gene.Symbol"] == "Agt"),
               which(GN379[,"Gene.Symbol"] == "Agtr1a"),
               which(GN379[,"Gene.Symbol"] == "Agtr1b")), c("Gene.Symbol")]
# ESR1 <- Estrogen receptor 1
ESR <- GN379[c(which(GN379[,"Gene.Symbol"] == "Esr1"),
               which(GN379[,"Gene.Symbol"] == "Esr2"),
               which(GN379[,"Gene.Symbol"] == "Esrrg")), c("Gene.Symbol")]
# APP <- Amyloid precursor protein (Alzheimer Disease)
APP <- GN379[c(which(GN379[,"Gene.Symbol"] == "App"),
               which(GN379[,"Gene.Symbol"] == "Saa1"),
               which(GN379[,"Gene.Symbol"] == "Saa1"),
               which(GN379[,"Gene.Symbol"] == "Apcs")), c("Gene.Symbol")]

# LEP <- Leptin (Obesity)
LEP <- GN379[c(which(GN379[,"Gene.Symbol"] == "Lep"),
               which(GN379[,"Gene.Symbol"] == "Lepr")), c("Gene.Symbol")]
# INS <- Insulin (Diabetes mellitus)
INS <- GN379[c(which(grepl("Insr", GN379[,"Gene.Symbol"])), 
               which(grepl("Ide", GN379[,"Gene.Symbol"])),
               which(grepl("Ins1", GN379[,"Gene.Symbol"])),
               which(grepl("Ins2", GN379[,"Gene.Symbol"])),
               which(grepl("Irs1", GN379[,"Gene.Symbol"])),
               which(grepl("Irs2", GN379[,"Gene.Symbol"]))), c("Gene.Symbol")]
# NPPB <- Brain Natriuretic Peptide (Cardiovascular disease)
NPP <- GN379[which(grepl("natriuretic", GN379[,"Description"])), c("Gene.Symbol")]

# Cell Cycle genes
CC <- GN379[c(which(grepl("cyclin", GN379[,"Description"]) & 
             !grepl("-like", GN379[,"Description"]) &
             !grepl("-related", GN379[,"Description"]) &
             !grepl("associated", GN379[,"Description"]) &
             !grepl("recycling", GN379[,"Description"]) &
             !grepl("-dependent", GN379[,"Description"]) )), c("Gene.Symbol")]

#Tumor
TUM <- GN379[c(which(grepl("tumor", GN379[,"Description"]) & 
             !grepl("-like", GN379[,"Description"]) &
             !grepl("superfamily", GN379[,"Description"]) &
             !grepl("associated", GN379[,"Description"]) &
             !grepl("open reading", GN379[,"Description"]) &
             !grepl("-dependent", GN379[,"Description"]) )), c("Gene.Symbol")]
             
highImpact <- unique(c(CD, IL, H2, ACE, ESR, APP, LEP, INS, NPP, CC, TUM))

phenotypes <- GN379[which(GN379[,"Gene.Symbol"] %in% highImpact),]
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
if(!file.exists("QTLs.txt")) {
  QTLs <- matrix(NA, nrow(phenotypes), nrow(genotypes))
  colnames(QTLs) <- rownames(genotypes)

  for(x in rownames(genotypes)){
    aovModel <- aov(apply(phenotypes, 1, as.numeric) ~ as.character(genotypes[x,]))
    QTLs[,x] <- unlist(lapply(summary(aovModel),function(x){ return(x[[5]][1]) }))
  }

  write.table(QTLs, file = "QTLs.txt", sep="\t", quote=FALSE)
}else{
  QTLs <- read.table("QTLs.txt")
}

# Write out the data as a cross object for R/qtl
phenotypes <- rbind(phenotypes, sex=rep("M", ncol(phenotypes)))
write.table(cbind(NA,NA, phenotypes), "cross.csvr", sep=",", quote=FALSE, col.names=FALSE, row.names=TRUE, na="")
write.table(cbind(map[,1],map[,2], genotypes), "cross.csvr", sep=",", append=TRUE, quote=FALSE, col.names=FALSE, row.names=TRUE, na="")
library(qtl)
cross <- read.cross("csvr", file="cross.csvr")

# Write out the data as a cross object for R/qtl
library(ctl)

cross <- fill.geno(cross)

if(!file.exists("CTLs_p.txt")) {
  CTLs <- matrix(NA, nrow(phenotypes), nrow(genotypes))
  colnames(CTLs) <- rownames(genotypes)

  cat("", file="CTLs_p.txt")
  cat("", file="CTLs_int.txt")
  for(x in 315:(nphe(cross)-1)){
    res   <- CTLscan.cross(cross, phenocol = x)
    scores <- apply(res[[1]]$ctl, 1, max)
    significant <- which(apply(res[[1]]$ctl, 2, max) > 2)                                     # Other phenotypes causing the CTL
    for(y in significant){
      gphes <- names(which(res[[1]]$ctl[,y] > 2))                                             # Markers causing the CTL
      locs  <- unique(cbind(map[gphes,1], round(as.numeric(map[gphes,3])/5) * 5))             # Approx location
      interacts <- cbind(annotation[x,1],annotation[y,1],locs, mean(res[[1]]$ctl[gphes,y]), "\n")
      cat(t(interacts), file="CTLs_int.txt", append = TRUE)
    }
    cat(paste0(paste0(c(x, scores),collapse="\t"), "\n"), file="CTLs_p.txt", append = TRUE)
    CTLs[x,] <- scores
    if(x %% 4 == 0) cat(x, "\n")
  }
  write.table(CTLs, file = "CTLs.txt", sep="\t", quote=FALSE)
}else{
  CTLs <- read.table("CTLs_p.txt", row.names=1)
}

haveQTL <- which(apply(-log10(QTLs), 1, max) > 5)
haveCTL <- which(apply(CTLs, 1, max) > 5)

dev.off()
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
  cat(annotation[x,1], max(-log10(QTLs[x,])), "\n")
}

for(x in haveCTL) {
  cat(annotation[x,1], max(CTLs[x,]), "\n")
}
