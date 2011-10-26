require(QCL)
setwd("E:/GBIC/YEAST")
probeannot <- read.csv("probeAnnotation.txt",sep="\t",header=TRUE,row.names=1)
cross <- read.cross("csvr",file="yeast_brem_cross.csv",geno=c(0,1))
cross <- convert2riself(cross)

phenotypes <- apply(pull.pheno(cross),2,as.numeric)
genotypes <- pull.geno(cross)

QCL <- QCL.scan(cross, 1:3, verbose=T)

setwd("e:/")
nodes <- QCLscanToSIF(QCL,0.7)
write.nodeAttributeFile(nodes,probeannot)
genes <- nodesToGenes(nodes,probeannot)
pathwaydata <- download.KEGG.annotation(genes,verb=T)