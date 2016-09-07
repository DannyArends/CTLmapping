#
# Analysis of GN111 BxD data, and selection of genes used in the analysis
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

GN111 <- read.csv(gzfile("GN111/GN111_MeanDataAnnotated_rev081815.txt.gz"), skip = 33, header = TRUE, sep="\t", colClasses="character")
GN111[1:10,]

GN111 <- GN111[- which(grepl("AFFX", GN111[,"Gene.Symbol"])), ]
GN111 <- GN111[- which(grepl("Rik", GN111[,"Gene.Symbol"])), ]
GN111 <- GN111[- which(grepl("Gm", GN111[,"Gene.Symbol"])), ]

# CD <- CD4/CD8 antigen // 126 genes
CD <- unique(GN111[c(which(grepl("Cd", GN111[,"Gene.Symbol"]) & 
              grepl("antigen", GN111[,"Description"]) & 
             !grepl("-like", GN111[,"Description"]) & 
             !grepl("cyclin", GN111[,"Description"]) & 
             !grepl("RIKEN", GN111[,"Description"]) & 
             !grepl("cadherin", GN111[,"Description"]) & 
             !grepl("homolog", GN111[,"Description"]) & 
             !grepl("division", GN111[,"Description"]) & 
             !grepl("family", GN111[,"Description"]))), c("Gene.Symbol")])
# IL <- Interleukines
IL <- unique(GN111[c(which(grepl("interleukin ", GN111[,"Description"]) & 
             !grepl("-like", GN111[,"Description"]) & 
             !grepl("-containing", GN111[,"Description"]) & 
             !grepl("binding protein", GN111[,"Description"]) & 
             !grepl("accessory", GN111[,"Description"]) & 
             !grepl("antagonist", GN111[,"Description"]) & 
             !grepl("transducer", GN111[,"Description"]) & 
             !grepl(", alpha", GN111[,"Description"]) & 
             !grepl(", beta", GN111[,"Description"]) & 
             !grepl("enhancer", GN111[,"Description"]) & 
             !grepl("family", GN111[,"Description"]))), c("Gene.Symbol")])
# HLA <- Major histocompatibility complex (MHC / H2)
H2 <- unique(GN111[c(which(grepl("H2-", GN111[,"Gene.Symbol"]))), c("Gene.Symbol")])
# ACE <-  Angiotensin-converting enzyme
ACE <- unique(GN111[c(which(GN111[,"Gene.Symbol"] == "Ace"),
               which(GN111[,"Gene.Symbol"] == "Ace2"),
               which(GN111[,"Gene.Symbol"] == "Agt"),
               which(GN111[,"Gene.Symbol"] == "Zbtb16"),
               which(GN111[,"Gene.Symbol"] == "Agtr1a"),
               which(GN111[,"Gene.Symbol"] == "Agtr1b")), c("Gene.Symbol")])
# ESR1 <- Estrogen receptor 1
ESR <- unique(GN111[c(which(GN111[,"Gene.Symbol"] == "Esr1"),
               which(GN111[,"Gene.Symbol"] == "Esr2"),
               which(GN111[,"Gene.Symbol"] == "Gh"),
               which(grepl("Igf", GN111[,"Gene.Symbol"])),
               which(GN111[,"Gene.Symbol"] == "Igf1"),
               which(GN111[,"Gene.Symbol"] == "Esrrg")), c("Gene.Symbol")])
# APP <- Amyloid precursor protein (Alzheimer Disease)
APP <- unique(GN111[c(which(GN111[,"Gene.Symbol"] == "App"),
               which(GN111[,"Gene.Symbol"] == "Saa1"),
               which(GN111[,"Gene.Symbol"] == "Saa1"),
               which(GN111[,"Gene.Symbol"] == "Apcs")), c("Gene.Symbol")])

# LEP <- Leptin (Obesity) / INS <- Insulin (Diabetes mellitus)
OBE <- unique(GN111[c(which(grepl("Insr", GN111[,"Gene.Symbol"])), 
               which(grepl("Ide", GN111[,"Gene.Symbol"])),
               which(grepl("Ins1", GN111[,"Gene.Symbol"])),
               which(GN111[,"Gene.Symbol"] == "Lep"),
               which(GN111[,"Gene.Symbol"] == "Lepr"),
               which(grepl("Ins2", GN111[,"Gene.Symbol"])),
               which(grepl("Insig2", GN111[,"Gene.Symbol"])),
               which(grepl("Negr1", GN111[,"Gene.Symbol"])),
               which(grepl("Fto", GN111[,"Gene.Symbol"])),
               which(grepl("Bdnf", GN111[,"Gene.Symbol"])),
               which(grepl("Pcsk1", GN111[,"Gene.Symbol"])),
               which(grepl("Irs1", GN111[,"Gene.Symbol"])),
               which(grepl("Npc1", GN111[,"Gene.Symbol"])),
               which(grepl("Bbs7", GN111[,"Gene.Symbol"])),
               which(grepl("Ccna2", GN111[,"Gene.Symbol"])),
               which(grepl("Irs2", GN111[,"Gene.Symbol"]))), c("Gene.Symbol")])
# NPPB <- Brain Natriuretic Peptide (Cardiovascular disease)
NPP <- unique(GN111[c(which(grepl("natriuretic", GN111[,"Description"])),
                      which(grepl("Comt", GN111[,"Gene.Symbol"])),
                      which(grepl("Npr", GN111[,"Gene.Symbol"]))), c("Gene.Symbol")])

# Cell Cycle genes
CC <- unique(GN111[c(which(grepl("cyclin", GN111[,"Description"]) & 
             !grepl("-like", GN111[,"Description"]) &
             !grepl("-related", GN111[,"Description"]) &
             !grepl("associated", GN111[,"Description"]) &
             !grepl("recycling", GN111[,"Description"]) &
             !grepl("-dependent", GN111[,"Description"]) )), c("Gene.Symbol")])

#Tumor
TUM <- unique(GN111[c(which(grepl("tumor", GN111[,"Description"]) & 
             !grepl("-like", GN111[,"Description"]) &
             !grepl("LOC", GN111[,"Gene.Symbol"]) &
             !grepl("superfamily", GN111[,"Description"]) &
             !grepl("associated", GN111[,"Description"]) &
             !grepl("open reading", GN111[,"Description"]) &
             !grepl("RIKEN", GN111[,"Description"]) &
             !grepl("EST", GN111[,"Description"]) &
             !grepl("-dependent", GN111[,"Description"]) )), c("Gene.Symbol")])
             
#HSP
HSP <- GN111[c(which(grepl("heat shock", GN111[,"Description"]) & 
             !grepl("-like", GN111[,"Description"]) )), c("Gene.Symbol")]

# Which group does a gene belong to
whichGroup <- function(x){
  if(x %in% CD) return("CD")
  if(x %in% IL) return("Il")
  if(x %in% H2) return("H2")
  if(x %in% ACE) return("ACE")
  if(x %in% ESR) return("ESR")
  if(x %in% APP) return("APP")
  if(x %in% OBE) return("OBE")
  if(x %in% NPP) return("NPP")
  if(x %in% CC) return("CC")
  if(x %in% TUM) return("TUM")
  if(x %in% HSP) return("HSP")
}

highImpact <- unique(c(CD, IL, H2, ACE, ESR, APP, OBE, NPP, CC, TUM, HSP))
cat("",file = "genes.txt")
for(x in highImpact){
  cat(x, "\t", whichGroup(x), "\n",file = "genes.txt", sep="", append=TRUE)
}
