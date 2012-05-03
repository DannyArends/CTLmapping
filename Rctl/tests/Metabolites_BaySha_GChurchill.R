#Analysis script for the Arabidopsis Bay x Sha Matabolites
# Paper senior author: Gary Churchill 2012, Published in: Plos
#
# (C) 2012 Danny Arends - GBIC - University of Groningen
#
# Start by either: Setwd("C:/github/CTLmapping/Rctl/tests")
# Files are distributed with the Rctl package in Rctl/tests

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

#Helper functions: 
# - Get a CTLobject's name
# - Remove the diagonal from a matrix
# - Chromosome edge locations
# - Color range for plots (Red-Black-Blue)
name <- function(CTLobject){ return(attr(CTLobject$ctl,"name")); }
remove.diag <- function(x){ return(x*lower.tri(x) + x*upper.tri(x)); }
get.chredges <- unlist(lapply(unique(map_info[,1]),function(x){max(which(map_info[,1]==x));}))
up <- abs(seq(-2,-0,0.1))/2 ; dw <- seq(0.1,2,0.1)/2
redblue <- c(rgb(up,0,0), rgb(0,0,dw))

#Load the library and scan the data
library(ctl)
ctls <- CTLscan(genotypes, metabolites, pheno.col=1:10, genotype.values=c("A","B"), n.perm = 25)

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

### Section: Input data quality control ###

png("Genetic_map.png",width=2000,height=1000)
  image(1:ncol(genotypes),1:ncol(genotypes),cor(genotypes=="A"),
        main="Genetic map",xlab="Markers",ylab="Markers",col=gray.colors(75))
  abline(h=(get.chredges+.5),col="white")
  abline(v=(get.chredges+.5),col="white"); box()
dev.off()

#Create the Two correlation matrices
cor_metabolites <- cor(metabolites,use='pair')
rownames(cor_metabolites) <- colnames(metabolites)
colnames(cor_metabolites) <- colnames(metabolites)
cor_individuals <- cor(t(metabolites),use='pair')
rownames(cor_individuals) <- rownames(metabolites)
colnames(cor_individuals) <- rownames(metabolites)

#Create a range of colors

png("Cor_Metabolites.png",width=2000,height=1000)
  image(x=1:ncol(metabolites),y=1:ncol(metabolites),cor_metabolites ,breaks=seq(-1,1.05,0.05),col=redblue)
  abline(h=(1:ncol(metabolites)+.5),col="white",lwd=0.1)
  abline(v=(1:ncol(metabolites)+.5),col="white",lwd=0.1); box()  
dev.off()

png("Cor_Individuals.png",width=1000,height=600)
  image(x=1:nrow(metabolites),y=1:nrow(metabolites),cor_individuals ,breaks=seq(-1,1.05,0.05),col=redblue)
dev.off()

png("Hist_CorIndividuals.png",width=1000,height=600)
  hist(cor_metabolites, breaks=seq(-1,1,0.1), col=rgb(1,0,0,0.5),freq=FALSE)
  hist(cor_individuals ,breaks=seq(-1,1,0.1), add=TRUE,col=rgb(0,0,1,0.5),freq=FALSE)
dev.off()

#Helper function to get the top-correlated metabolites
top.correlated <- function(x){
  ret <- t(apply(remove.diag(x),1,function(r){
    id <- which.max(abs(r))
    return(c(names(r)[id],id,r[id]))
  }))
  colnames(ret) <- c("top.correlated","id","correlation")
  return(ret)
}

top.correlated(cor_metabolites)[1:10,]
