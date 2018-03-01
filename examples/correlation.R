#
# R code to generate input and output files to check matrix*matrix pearson correlation
#

set.seed(101)             # Fix the seed for reproducable random number generation
n = 200                   # n = number of individuals, it's the shared dimension between A and B
nMeasurementsA = 100      # Number of measurements in A, is going to be the first dimension of C
nMeasurementsB = 50       # Number of measurements in B, is going to be the second dimension of C

A <- round(matrix(runif(n * nMeasurementsA), n, nMeasurementsA),3)   # 100 measurements, 200 individuals rounded down to 3 sigificant digits
B <- round(matrix(runif(n * nMeasurementsB), n, nMeasurementsB),3)   #  50 measurements, 200 individuals rounded down to 3 sigificant digits

A[sample(n * nMeasurementsA, 50)] <- 3   # Introduce missing data
B[sample(n * nMeasurementsB, 50)] <- 3   # Introduce missing data

# Write the input A for the C code
write.table(t(A), file = "inputA.txt", sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

# Write the input B for the C code
write.table(t(B), file = "inputB.txt", sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)

A[A == 3] <- NA  # Change the coding of the missing data in A to NA
B[B == 3] <- NA  # Change the coding of the missing data in B to NA

# calculate correlation for each column in A the correlation with matrix B
C <- apply(A, 2, function(columnA, mB) {
  return(cor(columnA, mB, method="pearson", use = "pair"))
}, B)

# Write the output C matrix, to compare the results of the C code against
write.table(round(t(C),12), file = "outputC.txt", sep = "\t", quote = FALSE, row.names=FALSE, col.names=FALSE)
q("no")

