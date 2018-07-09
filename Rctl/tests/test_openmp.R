#
# Testing reference implementation (R) versus unrolled loop (c) and openmp version (c)
# copyright (c) 2016-2020 - Danny Arends, Pjotr Prins, Gudrun Brockmann and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

library(ctl)

n.perm <- 20
n.ind <- 200
n.phe <- 200

time.ref <- 0
time.cor <- 0
time.omp <- 0

time.all <- matrix(NA, 0, 3, dimnames = list(c(), c("ref", "cor", "omp")))
for(i in 1:n.perm) {
  x <- runif(n.ind)                                   # Random x vector
  Y <- matrix(runif(n.ind * n.phe), n.ind, n.phe)     # Random Y matrix

  s.ref <- proc.time()[3]                             # Start time
  res.ref <- round(cor(x, Y), 6)                      # Run analysis
  time.ref <- time.ref + (proc.time()[3] - s.ref)     # Add time information

  s.cor <- proc.time()[3]
  res.cor <- round(correlation(x, Y), 6)
  time.cor <- time.cor + (proc.time()[3] - s.cor)

  s.omp <- proc.time()[3]
  res.omp  <- round(openmp(2, x, Y)$res, 6)
  time.omp <- time.omp + (proc.time()[3] - s.omp)

  cat(time.ref, time.cor, time.omp, "\n")

  time.all <- rbind(time.all, c(time.ref, time.cor, time.omp))

  if(!all(res.ref == res.cor)) stop("ref != cor\n")   # When results differ, stop
  if(!all(res.ref == res.omp)) stop("ref != mp\n")     # When results differ, stop
}

boxplot(time.all)
plot(time.all[,"ref"], time.all[,"cor"])
plot(time.all[,"ref"], time.all[,"omp"])
