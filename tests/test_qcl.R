library(qcl)

if(has_rqtl()){
  require(qtl)
  probeannot <- read.csv("small_annot.txt",sep="\t",header=TRUE,row.names=1)
  cross <- read.cross("csvr",file="brem_cross.csvr",geno=c("AA","AB"))
  cross <- convert2riself(cross)

  QCL <- QCLscan.cross(cross,1, verbose=T)
  image(QCL)
  plot(QCL)

  nodes <- QCLnetwork(QCL,0.55)

  write.nodeAttributeFile(nodes,probeannot)
  genes <- nodesToGenes(nodes,probeannot)
  pathwaydata <- download.KEGG.annotation(genes,verb=T)
}else{
  cat("Test skipped because of missing R/qtl package")
}