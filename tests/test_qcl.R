require(QCL)

probeannot <- read.csv("small_annot.txt",sep="\t",header=TRUE,row.names=1)
cross <- read.cross("csvr",file="brem_cross.csvr",geno=c("AA","AB"))
cross <- convert2riself(cross)

QCL <- QCL.scan(cross,1:3, verbose=T)

nodes <- QCLscanToSIF(QCL,0.55)

write.nodeAttributeFile(nodes,probeannot)
genes <- nodesToGenes(nodes,probeannot)
pathwaydata <- download.KEGG.annotation(genes,verb=T)
