#
# qtl.scan.R
#
# copyright (c) 2016-2020 - GBIC, Danny Arends, Pjotr Prins, Gudrun Brockmann, Rob Williams, Bruno Tesson and Ritsert C. Jansen
# last modified Nov, 2016
# first written Jan, 2011
# 
# R functions to do basic QTL mapping
#

# Maps a single trait QTL profile using a relatively 'slow' approach
QTLmapping <- function(genotypes, phenotypes, phenocol = 1, verbose = TRUE) {
  if(missing(genotypes) || is.null(genotypes)) stop("argument 'genotypes' is missing, with no default")
  if(missing(phenotypes)|| is.null(phenotypes)) stop("argument 'phenotypes' is missing, with no default")

  st <- proc.time()
  qtl <- rep(0, ncol(genotypes))
  tryCatch(
    qtl <- apply(genotypes, 2, function(x) {
                 -log10(anova(lm(phenotypes[,phenocol] ~ x))[[5]][1])
                }), error = function(e) e)
  if(verbose) cat("QTL mapping:", (proc.time()-st)[3], "seconds\n")
  invisible(qtl)
}

# end of qtl.scan.R

