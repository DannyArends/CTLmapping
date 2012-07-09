setwd("E:\\GBIC\\Lude CTL")
memory.limit(2000)

harmjanToLude <- function(hjnum = 32268){ koppel[which(koppel[,1]==hjnum),2] }

expdata    <- read.csv("out_RAW_eQTL.txt", sep="\t", row.names=1,check.names = FALSE)
transQTLs  <- read.csv("input/eQTLsFDR0.05.txt", header=TRUE, sep="\t",row.names=NULL)
koppel     <- read.csv("input/koppeltabel1.txt", header=TRUE, sep="\t",row.names=NULL)

probe_annot <- transQTLs[,c(5,17)]
probe_annot[,1] <- unlist(lapply(probe_annot[,1], harmjanToLude))

getGeneId <- function(name = "1400689"){
  as.character(unique(probe_annot[which(probe_annot[,1] == name),2]))
}

hjnames    <- unique(as.character(transQTLs[,5]))
allnames   <- unique(as.character(unlist(lapply(hjnames,harmjanToLude))))

snp_data <- read.table('genotypes/output.ped',header=FALSE,sep=" ",row.names=2,as.is=TRUE)[,-c(1:5)]
snp_descr <- read.table('genotypes/output.map',header=FALSE,sep=" ",row.names=2)

#Match to the genotypes
expdata <-  expdata[, match(rownames(snp_data), colnames(expdata))]
#Match to the transQTLs
expdata <-  expdata[match(allnames,rownames(expdata)),]

firstbase     <- snp_data[, seq(1,ncol(snp_data),2)]
secondbase    <- snp_data[, seq(2,ncol(snp_data),2)]

getSNP <- function(firstbase){
  if(levels(firstbase)[1]=='0') return(levels(firstbase)[-1])
  return(levels(firstbase))
}

createGenotype_N3 <- function(firstbase,secondbase, x){
  genotype   <- NULL
  snps       <- getSNP(as.factor(firstbase[,x]))
  mallele    <- names(which.max(table(unlist(firstbase[,1],secondbase[,1]))))
  for(b in 1:length(firstbase[,x])){
    if(as.character(firstbase[b,x])==as.character(secondbase[b,x])){
      if(as.character(firstbase[b,x])==snps[1]){
        genotype <- c(genotype,1)
      }else{
        genotype <- c(genotype,3)
      }
    }else{
      genotype <- c(genotype,2)
    }
  }
  return(genotype)
}

qtlmatrix <- NULL
for(snp in 1:ncol(firstbase)){
  cat('SNP', snp,'\n')
  genotype   <- createGenotype_N3(firstbase, secondbase, snp)
  models     <- aov(apply(expdata,1,as.numeric) ~ as.numeric(as.factor(genotype)))
  modelinfo  <- summary(models)
  qtls       <- -log10(as.numeric(unlist(lapply(modelinfo,"[",1,5))))
  qtlmatrix  <- cbind(qtlmatrix,qtls)
}

sign_snps   <- which(apply(qtlmatrix,2,max,na.rm=T) > 5)
sign_snps   <- sign_snps[-c(1:2)]
sign_probes <- which(apply(qtlmatrix,1,max,na.rm=T) > 5)

ctlmatrix <-   CTLprofiles(ctls,signi=2)
#image(1:153, 1:103, t(qtlmatrix[sign_probes,sign_snps]),col=c("white","gray","black"),breaks=c(0,5,10,1000),ylab='probe',xlab='snp',xaxt='n')
image(1:153, 1:103, t(ctlmatrix),col=c("white",gray.colors(2)[2:1],"black"),breaks=c(0,2,5,10,1000),ylab='probe',xlab='snp',xaxt='n')
box()
#chr_start <- c()
#chr_ends <- c()
#for(x in unique(snp_descr[sign_snps,1])){
#  chr_start <- c(chr_start,min(which(snp_descr[sign_snps,1]==x)))
#  chr_ends  <- c(chr_ends, max(which(snp_descr[sign_snps,1]==x)))
#}
abline(v=.5 + chr_ends)
for(x in 1:length(chr_start)){
  axis(1,at=(chr_start[x]+chr_ends[x])/2,paste(unique(aa[,2])[x]),cex=0.3)
}
abline(h=(1:200)+.5,col='white')
box()


plot.lodprofile <- function(aa, sign_probes, sign_snps, probe){
  plot(aa[probe,],t='h',xaxt='n',lwd=3)
  for(x in 1:length(chr_start)){
    axis(1,at=(chr_start[x]+chr_ends[x])/2,paste(unique(snp_descr[sign_snps,1])[x]),cex=0.3)
  }
  abline(v=.5 + chr_ends,col='red',lty=3)
}

plot.lodprofile <- function(qtlmatrix, sign_probes, sign_snps, probe){
  plot(qtlmatrix[sign_probes[probe],sign_snps],t='h',xaxt='n',lwd=3)
  for(x in 1:length(chr_start)){
    axis(1,at=(chr_start[x]+chr_ends[x])/2,paste(unique(snp_descr[sign_snps,1])[x]),cex=0.3)
  }
  abline(v=.5 + chr_ends,col='red',lty=3)
}

createGenotype_N2 <- function(firstbase,secondbase, x){
  genotype   <- NULL
  snps       <- getSNP(as.factor(firstbase[,x]))
  mallele    <- names(which.max(table(unlist(firstbase[,x],secondbase[,x]))))
  for(b in 1:length(firstbase[,x])){
    if(as.character(firstbase[b,x])==as.character(secondbase[b,x])){
      if(as.character(firstbase[b,x])==mallele){
        genotype <- c(genotype,1)
      }else{
        genotype <- c(genotype,NA)
      }
    }else{
      genotype <- c(genotype,2)
    }
  }
  cat(mallele,' ',length(which(genotype==1)),' ',length(which(genotype==2)),'\n')
  return(genotype)
}

genotypes_n3 <- NULL
genotypes_n2 <- NULL
for(snp in sign_snps){
  cat('SNP', snp,'\n')
  genotypes_n3   <- rbind(genotypes_n3, createGenotype_N3(firstbase, secondbase, snp))
  genotypes_n2   <- rbind(genotypes_n2, createGenotype_N2(firstbase, secondbase, snp))
}

write.table(file="hjg.txt",cbind(snp_descr[sign_snps,c(1,3)],genotypes_n2[,]),sep="\t",na="",col.names=FALSE)
write.table(file="hjpRAW.txt",cbind(NA,NA,expdata[sign_probes,]),sep="\t",na="",col.names=FALSE,quote=FALSE)

 library(ctl)
 ctls <- load.ctl("hjg.txt", "hjp.txt", "hjo")
 plot(ctls)


#Plot QTL and CTL
plot(aa[21,],t='h',col='blue',lwd=3,main="T21")
points(qtlmatrix[sign_probes[21],sign_snps], t='h',ylab='# QTL',lwd=2)

#Get all the probes that show a QTL on marker 22
rownames(expdata)[sign_probes[which(qtlmatrix[sign_probes,sign_snps[22]] >5)]]

#Get all significant CTLs
for(x in 1:length(resHJ)){
  mysign <- names(which(apply(resHJ[[x]]$l,1,function(x){any(x > 2)})))
  if(length(mysign) > 0) cat(x,' ',getGeneId(attr(resHJ[[x]]$ctl,'name')), '->', paste(lapply(mysign,getGeneId), collapse=', '),'\n')
}


qtls <- c("1820386", "7200044", "1400689", "7570196", "1470685", "1050292", "1770170", "6290747", "3850630",
  "1580025", "5560280", "130433",  "6040259", "5310437", "4880600", "5050086", "50136",   "1300431", 
  "630470",  "1050008", "2340577", "7550343", "3830327", "5810685")
unlist(lapply(qtls,getGeneId))

