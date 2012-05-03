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

ctl <- CTLscan(genotypes, metabolites, genotype.values=c("A","B"), n.perm = 100)

#Create comparison QTL / CTL heatmaps
png("Comparison_QTL_CTL.png",width=2000,height=1000)
  op <- par(mfrow=c(1,2))
  QTLimage(ctl)
  image(ctl)
dev.off()

