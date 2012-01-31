#plots for the paper
library(qcl)
qcl_result <- QCLscan(ath.metab$genotypes, ath.metab$phenotypes, n.perm=200)
setwd("E:\\github\\Rpackages\\QCLmapping\\inst\\methodPaper\\img")


pdf("ctl_qtl.pdf",width=10,height=7)
  plot.QCLscan2(qcl_result[[3]])
  text(54,10,"P1")
  arrows(54,9,y1=7,col="red",lwd=2)
  text(73,10,"P2")
  arrows(73,9,y1=7,col="red",lwd=2)
  text(93,10,"P3")
  arrows(93,9,97,y1=7,col="red",lwd=2)
dev.off()

pdf("ctl_pxp.pdf",width=10,height=7)
  op <- par(mar=c(7, 8, 4, 2) + 0.1)
  image(qcl_result,against="phenotypes")
dev.off()

QCLnetwork(qcl_result)

#Open cytoscape and create the last plot from network.sif