setwd("E:\\GBIC\\Lude CTL")
source('functions.r')
library(ctl)
memory.limit(2000)
transQTLs      <- read.csv("input/eQTLsFDR0.05.txt", header=TRUE, sep="\t",row.names=NULL)
koppel         <- read.csv("input/koppeltabel1.txt", header=TRUE, sep="\t",row.names=NULL)

#Prepare smaller subsets using the eQTLsFDR0.05
expdata_pca    <- read.csv("out_40PCA_noG.txt", sep="\t", row.names=1)
expdata_raw    <- read.csv("out_RAW.txt", sep="\t", row.names=1)


hjnames    <- as.character(transQTLs[,5])
allnames   <- unique(as.character(unlist(lapply(hjnames,harmjanToLude))))
expids_pca <- which(rownames(expdata_pca) %in% allnames)
expids_raw <- which(rownames(expdata_raw) %in% allnames)

expdata_pca  <- expdata_pca[expids_pca, ]
expdata_raw  <- expdata_raw[expids_raw, ]
namesexp_pca <- gsub(".","-",colnames(expdata_pca),fixed=T)
namesexp_raw <- gsub(".","-",colnames(expdata_raw),fixed=T)
namesexp_pca <- gsub("HT12_","-",namesexp_pca)
namesexp_raw <- gsub("HT12_","-",namesexp_raw)
colnames(expdata_pca) <- namesexp_pca
colnames(expdata_raw) <- namesexp_raw

write.table(expdata_pca,file='out_40PCA_noG_small.txt',quote=FALSE, sep='\t')
write.table(expdata_raw,file='out_RAW_small.txt'      ,quote=FALSE, sep='\t')

#Load back the small dataset and check row and colnames matches
in_pca  <- read.table("out_40PCA_noG_small.txt",check.names = FALSE)
in_raw  <- read.table("out_RAW_small.txt",check.names = FALSE)

raw_ids <- match(colnames(in_pca),colnames(in_raw))
in_raw  <- in_raw[,raw_ids]

in_pca  <- in_pca[rownames(in_raw),] # we cannot match 1 in the raw that is in PCA ?

#Map QTLs to select a smaller subset of genes
snp_data  <- read.table('genotypes/output.ped',header=FALSE,sep=" ",row.names=2,as.is=TRUE)[,-c(1:5)]
snp_descr <- read.table('genotypes/output.map',header=FALSE,sep=" ",row.names=2)

in_pca <-  in_pca[, match(rownames(snp_data), colnames(in_pca))]
in_raw <-  in_raw[, match(rownames(snp_data), colnames(in_raw))]

firstbase     <- snp_data[, seq(1,ncol(snp_data),2)]
secondbase    <- snp_data[, seq(2,ncol(snp_data),2)]

qtl_raw <- mapQTLs(in_raw, firstbase, secondbase)
colnames(qtl_raw) <- rownames(snp_descr)
rownames(qtl_raw) <- rownames(in_pca)
qtl_pca <- mapQTLs(in_pca, firstbase, secondbase)
colnames(qtl_pca) <- rownames(snp_descr)
rownames(qtl_pca) <- rownames(in_pca)

probe_s <- sort(unique(c(which(apply(qtl_raw,1,function(x){any(x>5,na.rm=T)})), which(apply(qtl_pca,1,function(x){any(x>5,na.rm=T)})))))
snps_s  <- sort(unique(c(which(apply(qtl_pca,2,function(x){any(x>5,na.rm=T)})), snps_s_raw <- which(apply(qtl_pca,2,function(x){any(x>5,na.rm=T)})))))

image(t(qtl_raw[probe_s,snps_s]),col=c('white','black'),breaks=c(0,5,1000))
image(t(qtl_pca[probe_s,snps_s]),col=c('white','black'),breaks=c(0,5,1000))

genotypes_n2 <- NULL
for(snp in 1:ncol(firstbase)){
  cat('SNP', snp,'\n')
  genotypes_n2   <- rbind(genotypes_n2, createGenotype_N2(firstbase, secondbase, snp))
}

colnames(genotypes_n2) <- rownames(firstbase)
rownames(genotypes_n2) <- rownames(snp_descr)

enough_1  <- names(which(apply(genotypes_n2,1,function(x){sum(x==1,na.rm=T)}) > 200))
enough_2  <- names(which(apply(genotypes_n2,1,function(x){sum(x==2,na.rm=T)}) > 200))
good_snps <- c(enough_1,enough_2)[which(duplicated(c(enough_1,enough_2)))]

genotypes <- genotypes_n2[good_snps,]

write.table(cbind(snp_descr[good_snps,c(1,3)],genotypes), file='new_geno.txt',col.names=FALSE,na="",sep="\t")
write.table(cbind(NA,NA,in_pca),file='pheno_pca.txt',col.names=FALSE,na="",sep="\t",quote=FALSE)
write.table(cbind(NA,NA,in_raw),file='pheno_raw.txt',col.names=FALSE,na="",sep="\t",quote=FALSE)

#We can restart from here after running mapctl.exe
setwd("E:\\GBIC\\Lude CTL")
source('functions.r')
memory.limit(2000)
library(ctl)
good_snps <- as.character(unlist(read.table("good_snps.txt")))
snp_descr <- read.table('genotypes/output.map',header=FALSE,sep=" ",row.names=2)

raw_ctls <- ctl.load("new_geno.txt", "pheno_raw.txt", "ctl_raw",verbose=TRUE)
pca_ctls <- ctl.load("new_geno.txt", "pheno_pca.txt", "ctl_pca",verbose=TRUE)
cor_ctls <- ctl.load("new_geno.txt", "pheno_cor.txt", "ctl_cor",verbose=TRUE,to=150)

#op <- par(mfrow = c(2,1))
#internal.image(qtl_raw[1:150,good_snps],breaks= c(0, 3, 5, 10, 15, 10000), snp_descr=snp_descr[good_snps,])
image(pca_ctls, significance=1.1, grid.col="white", snp_descr=snp_descr[good_snps,])


idx_raw <- which(apply(qtl_raw[,good_snps],1,max) > 3)
idx_pca <- which(apply(qtl_pca[,good_snps],1,max) > 3)

png("qtls_raw.png",w=1800,h=1200)
op <- par(cex= 2)
internal.image(qtl_raw[idx_raw,good_snps],breaks= c(0, 3, 5, 10, 15, 10000), snp_descr=snp_descr[good_snps,],main='QTLs on RAW data')
dev.off()

png("qtls_pca.png",w=1800,h=1200)
op <- par(cex= 2)
internal.image(qtl_pca[idx_pca,good_snps],breaks= c(0, 3, 5, 10, 15, 10000), snp_descr=snp_descr[good_snps,],main='QTLs on PCA corrected data')
dev.off()

profiles <- plot(ctls,sign=0.05, snp_descr=snp_descr[good_snps,])
for(x in rownames(profiles)){cat(which(ctl.names(ctls)==x)-1,"\n")}

png("ctls_raw.png",w=1800,h=1200)
op <- par(cex= 2)
image(raw_ctls, significance=0.05, grid.col="white", snp_descr=snp_descr[good_snps,])
dev.off()

png("ctls_pca.png",w=1800,h=1200)
op <- par(cex= 2)
image(pca_ctls, significance=0.05, grid.col="white", snp_descr=snp_descr[good_snps,])
dev.off()


names <- rownames(CTLprofiles(raw_ctls))

GeneNames <- transQTLs[,c(5,17)]
genename <- function(x){
  unlist(lapply(x, function(e){ as.character(unique(GeneNames[which(GeneNames[,1]==LudeToHJ(e)),2]))}))
}

for(x in cor_ctls){
  idx <- which(apply(x$l,1,max) > -log10(0.05))
  if(length(idx) > 0) cat(genename(ctl.name(x)),"->",genename(names(idx)), "(",idx,",",names(idx),")\n")
}



png("sum_ctls_pca.png",w=1800,h=1000)
op <- par(cex= 2)
  plot(apply(aa,2,sum),t='l',xaxt='n',ylab="Summed LOD", xlab="Marker")
  addChromosomeLines(snp_descr[good_snps,],col='blue')
dev.off()

png("sum_ctls_raw.png",w=1800,h=1000)
op <- par(cex= 2)  
  plot(apply(aa,2,sum),t='l',xaxt='n',ylab="Summed LOD", xlab="Marker")
  addChromosomeLines(snp_descr[good_snps,],col='blue')
dev.off()

names(which(aa[,"rs2540917"] > 1))

#COHORT data correction
cohort <- as.numeric(as.factor(substr(colnames(in_raw),0,3)))
plot(-log10(unlist(lapply(summary(aov(apply(in_raw,1,as.numeric) ~ cohort)),"[",1,5))),t='l')

correctionmatrix <- NULL
for(p in 1:nrow(in_raw)){
  difmeans <- mean(as.numeric(in_raw[p,])) - unlist(lapply(unique(cohort),function(x){mean(as.numeric(in_raw[p,cohort==x]))}))
  traitcorrection <- rep(0,ncol(in_raw))
  for(x in unique(cohort)){
    traitcorrection[cohort==x] <- difmeans[x]
  }
  correctionmatrix <- rbind(correctionmatrix,traitcorrection)
}
in_cor <- in_raw - correctionmatrix
write.table(cbind(NA,NA,in_cor),file='pheno_cor.txt',col.names=FALSE,na="",sep="\t",quote=FALSE)

plot(-log10(unlist(lapply(summary(aov(apply(in_raw-correctionmatrix,1,as.numeric) ~ cohort)),"[",1,5))),t='l')


cc <-NULL
for(x in 1:499){cc <- c(cc,cor(as.numeric(in_raw[x,]),as.numeric(in_cor[x,])))}

