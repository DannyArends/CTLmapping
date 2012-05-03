setwd("E:/GBIC/Metabolomics")
metabolites <- read.table("expressions.txt",sep="\t",row.names=1,header=TRUE)
genotypes   <- read.table("genotypes.txt",row.names=1,header=TRUE)
map_info <- apply(t(genotypes[1:3,]),2,as.numeric)
rownames(map_info) <- colnames(genotypes)
genotypes <- genotypes[-c(1:3),]
rownames(genotypes) <- as.numeric(rownames(genotypes))

ids <- match(rownames(metabolites),rownames(genotypes))
genotypes <- genotypes[ids,]

metab_means <- grep("Mean",colnames(metabolites))
metabolites <- metabolites[,metab_means]

metabolites[1:10,1:10]
genotypes[1:10,1:10]
map_info[1:10,1:3]

library(ctl)

ctl <- CTLscan(genotypes, metabolites, n.perm = 1500, method = "kendall",genotype.values=c("A","B"))

png("Comparison_QTL_CTL.png",width=2000,height=1000)
  op <- par(mfrow=c(1,2))
  QTLimage(ctl)
  image(ctl)
dev.off()
