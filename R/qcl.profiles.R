#
# qcl.profiles.R
#
# copyright (c) 2011 Danny Arends and Ritsert C. Jansen
# last modified Oct, 2011
# first written Oct, 2010
# 
#

#Create the 2 possible QCL matrices
#Phenotypes versus markers
#Phenotypes versus phenotypes
QCLprofiles <- function(QCLscan, qcl.threshold=0.4, against = c("markers","phenotypes")){
  mymatrix <- NULL
  phenotypenames <- NULL
  if(against[1] == "markers"){
    for(p in 1:length(QCLscan)){
      mymatrix <- rbind(mymatrix,apply(QCLscan[[p]],2,function(x){length(which(x > qcl.threshold))}))
      phenotypenames <- c(phenotypenames,attr(QCLscan[[p]],"name"))
    }
    rownames(mymatrix) <- phenotypenames
    colnames(mymatrix) <- colnames(QCLscan[[1]])
    return(mymatrix)
  }
  if(against[1] == "phenotypes"){
    targets <- NULL
    for(p in 1:length(QCLscan)){
      targets <- unique(c(targets,unique(as.character(unlist(apply(QCLscan[[p]],2,function(x){names(which(x > qcl.threshold))}))))))
      phenotypenames <- c(phenotypenames,attr(QCLscan[[p]],"name"))
    }
    mymatrix <- matrix(0,length(QCLscan),length(targets))
    colnames(mymatrix) <- targets
    for(p in 1:length(QCLscan)){
      current_table <- table(as.character(unlist(apply(QCLscan[[p]],2,function(x){names(which(x > qcl.threshold))}))))
      mymatrix[p,names(current_table)] <- current_table
    }
    rownames(mymatrix) <- phenotypenames
    colnames(mymatrix) <- targets
    mymatrix
  }
}
