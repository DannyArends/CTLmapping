setwd("/lustre/medusa/dannya")

phenotypes <- read.csv("GN320_MeanDataAnnotated_rev081815.txt", sep="\t", skip=32, header=TRUE, check.names=FALSE)
phenotypes <- phenotypes[, -441]
write.table(cbind("", "", "", phenotypes[,-c(1:14)]), file="GN302_pheno.csv", sep="\t", row.names=TRUE, col.names=FALSE,na="-")

setwd("/lustre/medusa/dannya")
map <- read.csv("HLC.map", sep="\t", header=FALSE, check.names=FALSE)
lines <- readLines("HLC.ped")

genotypes <- matrix(NA, 358520, length(lines))
for(ind in 1:length(lines)){
    l1 <- strsplit(lines[ind], "\t")[[1]]
    l2 <- l1[-1]
    geno <- rep(NA, length(l2)/2)
    cnt <- 1
    for(x in seq(1, length(l2)-1, 2)){
        geno[cnt] <- paste0(l2[x],l2[x+1], collapse="")
        cnt <- cnt + 1
    }
    geno[geno == "00"] <- NA
    genotypes[,ind] <- geno
    ind <- ind + 1
    cat("Done", ind, "\n")
}
write.table(genotypes, "GN302_geno.csv", sep="\t")

genotypes <- read.csv("GN302_geno.csv", sep="\t")

names <- c()
for(ind in 1:length(lines)){
    l1 <- strsplit(lines[ind], "\t")[[1]]
    names <- c(names, l1[1])
}

rownames(genotypes) <- map[,2]
mapchr1 <- map[which(map[,4] > 0 & map[,1] == 1), -3][,c(2,1,3)]
genochr1 <- genotypes[mapchr1[,1],]

write.table(genochr1, "GN302_geno_chr1.csv", sep="\t")
write.table(mapchr1, "GN302_map_chr1.csv", sep="\t")

genochr1num <- t(apply(genochr1,1,function(x){
    return(as.numeric(as.factor(as.character(x))))
}))
write.table(cbind(mapchr1[,c(2,3)], genochr1), file="GN301_genotypes.csv",sep="\t", row.names=TRUE, col.names=FALSE)

pedigree <- read.csv("HLC.ped", nrow = 4)
