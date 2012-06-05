setwd("e:/GBIC/MIT")
library(qtl)
library(ctl)
phenotypes <- read.csv("REF.GCT", sep="\t",header=T,row.names=1)
descriptions <- phenotypes[,1]
phenotypes <- phenotypes[,-1]
genodata <- read.csv("AXBfullmap1.csv",sep=",",header=T)
unique_genotypes <- genodata[!duplicated(genodata[,4:29]),]

pheno4cross <- cbind(rownames(phenotypes),NA,NA,phenotypes)
colnames(pheno4cross) <- colnames(unique_genotypes)
crosscsvr   <- rbind(pheno4cross,unique_genotypes)

write.table(crosscsvr,"AXBuniquemap.csv",sep=",",col.names=FALSE,row.names=FALSE,quote=FALSE,na="")

cross <- read.cross("csvr",file="AXBuniquemap.csv",genotypes=c("AA","H","BB"))
class(cross)[1] <- "riself"

#Minor fix
cross <- drop.markers(cross,"rs13478093")
filled_cross <- fill.geno(cross)
filled_cross <- calc.genoprob(filled_cross)

#We have a cross, now the analysis
if(file.exists("scanone_qtls.csv")){
  qtls <- read.csv("scanone_qtls.csv",row.names=1)
  qtls <- apply(qtls,2,as.numeric)
}else{
  qtls <- NULL
  for(x in seq(1,nphe(filled_cross),100)){
    if(x+99 > nphe(filled_cross)){
      m <- nphe(filled_cross)-2
    }else{ m <- x+99; }
    res <- scanone(filled_cross,pheno.col=x:m)
    qtls <- cbind(qtls,as.matrix(res[,3:ncol(res)]))
    cat("Done",x,"..",x+99,"\n")
  }
}

todo <- NULL
chrs <- rep(names(nmar(cross)),nmar(cross))
for(x in unique(chrs)){
  maxqtls     <- apply(qtls[chrs==x,],2,max)
  significant <- which(maxqtls > 5)
  cat("Chr",x,": ",length(significant),"\n")
  if(length(significant) > 1){
    image(1:nrow(qtls),1:length(significant),qtls[,significant],col=c("white",gray.colors(10)[10:1],"black"),main="QTL heatmap")
    box()
    todo <- c(todo,significant)
  }
}

filled_cross$pheno <- filled_cross$pheno[,todo]
qtls <- qtls[,todo]
ctls <- CTLscan.cross(filled_cross, pheno.col=1:250, have.qtl=t(qtls), n.perm=250)
