#Analysis script for the Arabidopsis Bay x Sha Matabolites
# Paper senior author: Gary Churchill 2012, Published in: Plos
#
# (C) 2012 Danny Arends - GBIC - University of Groningen
#   - Files are distributed with the Rctl package in Rctl/tests
# Start by: 
# setwd("C:/github/CTLmapping/Rctl/tests")

metabolites <- read.table("Metabolites_BaySha_GChurchill.txt",sep="\t",row.names=1,header=TRUE)
genotypes   <- read.table("Genotypes_BaySha.txt",row.names=1,header=TRUE)

#Create the genetic map from the first 3 lines og the genotypes file
map_info <- apply(t(genotypes[1:3,]),2,as.numeric)
rownames(map_info) <- colnames(genotypes)

#Genotype matrix by dropping the first 3 lines with map information
genotypes <- genotypes[-c(1:3),]
rownames(genotypes) <- as.numeric(rownames(genotypes))

#Match the genotype matrix to the phenotypes
ids <- match(rownames(metabolites),rownames(genotypes))
genotypes <- genotypes[ids,]

#Take only the mean expressions.
#How can we use the StdDev?
metab_means <- grep("Mean",colnames(metabolites))
metabolites <- metabolites[,metab_means]

#Print a piece of the data how does it look?
metabolites[1:5,1:10]
genotypes[1:5,1:10]
map_info[1:10,1:3]

#Load the library and scan the data
library(ctl)
source("Helper_Functions.r")
source("Basic_QC.R",local=TRUE,echo=TRUE)
ctls <- CTLscan(genotypes, metabolites, geno.enc=c("A","B"), pheno.col=1:10, n.perm = 25)

#Create comparison QTL / CTL heatmaps using QTLimage() and image.CTLscan()
png("Comparison_QTL_CTL.png",width=2000,height=1000)
  op <- par(mfrow=c(1,2)); op <- par(cex=1.8);
  QTLimage(ctls); image(ctls);
dev.off()

#Plot the individual CTL plot.CTLobject()
for(ctl in ctls){
  png(paste("CTL_",name(ctl),".png",sep=""),width=2000,height=1000); 
    op <- par(mfrow=c(1,2))
    op <- par(cex=1.8)
    plot(ctl)
    image(1:ncol(ctl$ctl),1:nrow(ctl$ctl),t(ctl$l),col=gray.colors(10)[10:1],ylab="Metabolites",xlab="Markers",main="Phenotype @ Marker x Phenotype matrix")
    abline(v=(get.chredges+.5),col="white");box()
  dev.off()
}
