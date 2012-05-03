#setwd("C:/github/CTLmapping/Rctl/tests")
metabolites <- read.table("Metabolites_BaySha_GChurchill.txt",sep="\t",row.names=1,header=TRUE)
genotypes   <- read.table("Genotypes_BaySha.txt",row.names=1,header=TRUE)
#Create A map from the first 3 lines
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

library(ctl)

ctls <- CTLscan(genotypes, metabolites, pheno.col=1:3, genotype.values=c("A","B"), n.perm = 5)

#Create comparison QTL / CTL heatmaps
png("Comparison_QTL_CTL.png",width=2000,height=1000)
  op <- par(mfrow=c(1,2)); QTLimage(ctls); image(ctls);
dev.off()

name <- function(ctl){ return(attr(ctl$ctl,"name")); }

#Plot the individual CTLs
for(ctl in ctls){
  png(paste("CTL_",name(ctl),".png",sep=""),width=2000,height=1000); plot(ctl); dev.off()
}

cor_metabolites <- cor(metabolites,use='pair')
rownames(cor_metabolites) <- colnames(metabolites)
colnames(cor_metabolites) <- colnames(metabolites)
cor_individuals <- cor(t(metabolites),use='pair')
rownames(cor_individuals) <- rownames(metabolites)
colnames(cor_individuals) <- rownames(metabolites)

redblue <- c(rgb(abs(seq(-1,0.11,0.1)),0,0), rgb(0,0,seq(0.11,1,0.1)))


png("Cor_Metabolites.png",width=2000,height=1000)
  image(x=1:ncol(metabolites),y=1:ncol(metabolites),cor_metabolites ,breaks=seq(-1,1,0.095),col=redblue)
dev.off()

png("Cor_Individuals.png",width=1000,height=600)
  image(x=1:nrow(metabolites),y=1:nrow(metabolites),cor_individuals ,breaks=seq(-1,1,0.095),col=redblue)
dev.off()

png("Hist_CorIndividuals.png",width=1000,height=600)
  hist(cor_metabolites, breaks=seq(-1,1,0.1), col=rgb(1,0,0,0.5),freq=FALSE)
  hist(cor_individuals ,breaks=seq(-1,1,0.1), add=TRUE,col=rgb(0,0,1,0.5),freq=FALSE)
dev.off()

remove.diag <- function(x){ return(x*lower.tri(x) + x*upper.tri(x)); }

top.correlated <- function(x){
  ret <- t(apply(remove.diag(x),1,function(r){
    id <- which.max(abs(r))
    return(c(names(r)[id],id,r[id]))
  }))
  colnames(ret) <- c("top.correlated","id","correlation")
  return(ret)
}

top.correlated(cor_metabolites)[1:10,]
