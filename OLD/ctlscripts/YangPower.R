library(ctl)

plot(c(0, 1.5),c(0,5),xlab="Dcor",ylab="LOD",t='n')
c_id <- 1
colorz <- c("red","orange","blue","green")

for(n.population in c(50)){
  n.phenotypes = 3

  VALS <- NULL
  for(x in 1:100){
    repeat{
      traits <- matrix(rnorm(n.population*n.phenotypes),n.population,n.phenotypes)
      marker <- round(runif(n.population),d=0)+1
      cat(".")
      if(dcor(t(t(marker)),traits,1,1,2)[3] > x/10000){
        cat("\n")
        break
      }
      if(dcor(t(t(marker)),traits,1,1,3)[3] > x/10000){
        cat("\n")
        break
      }
    }
    res <- CTLmapping(t(t(marker)), traits, n.perms = 10000, verbose = FALSE)
    VALS <- rbind(VALS, c(dcor(t(t(marker)),traits,1,1,2),res$ctl[2]))
    VALS <- rbind(VALS, c(dcor(t(t(marker)),traits,1,1,3),res$ctl[3]))
  }
  points(VALS[,c(3,4)], col=colorz[c_id], pch=19, cex=0.5)
  c_id <- c_id + 1
  cat("Done with population size:",n.population,"\n")
}

