#include "mapctl.h"

/* Perform a CTL scan and permutations on phenotype 'phenotype' */
double** mapctl(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                bool doperms, int nperms, bool verbose){

  info("Phenotype %d: Mapping", (phenotype+1));
  clvector* genoenc = getGenotypes(genotypes, false);
  size_t i;  
  double** ctls;
  double*  perms;
  double** scores = ctleffects(phenotypes, genotypes, phenotype, genoenc, 1, 1, verbose);
  if(!doperms){
    info(", toLOD\n");  // Exact calculation can be used
    ctls = toLODexact(scores, genoenc, genotypes.nmarkers, phenotypes.nphenotypes);
  }else{
    info(", Permutation");
    perms = permute(phenotypes, genotypes, phenotype, genoenc, 1, 1, nperms, false);
    info(", toLOD\n");
    ctls = toLOD(scores, perms, genotypes.nmarkers, phenotypes.nphenotypes, nperms);
    free(perms);
  }
  for(i = 0; i < genotypes.nmarkers; i++){ free(genoenc[i].data); }
  free(genoenc);
  freematrix((void**)scores, genotypes.nmarkers);
  return ctls;
}

/* Calculate the difference in correlation matrix for phenotype 'phenotype' */
double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                    clvector* genoenc, int alpha, int beta, bool verbose){

  size_t g, m, p,debug = 0;
  double** difcormatrix = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);

  for(m = 0; m < genotypes.nmarkers; m++){
    size_t ngenotypes = genoenc[m].nelements;

    clvector*  splits = (clvector*) calloc(ngenotypes, sizeof(clvector));
    double**   pheno1 = (double**)  calloc(ngenotypes, sizeof(double*));
    for(g = 0; g < ngenotypes; g++){
      splits[g] = which(genotypes.data[m], phenotypes.nindividuals, genoenc[m].data[g]);
      pheno1[g] = get(phenotypes.data[phenotype], splits[g]);
    }
    for(p = 0; p < phenotypes.nphenotypes; p++){
      if(p!=phenotype){
        double* cors  = newdvector(ngenotypes);
        int* nsamples = newivector(ngenotypes);
        for(g = 0; g < ngenotypes; g++){
          double* P2   = get(phenotypes.data[p], splits[g]);
          cors[g]      = correlation(pheno1[g], P2, splits[g].nelements, false);
          nsamples[g]  = splits[g].nelements;
          free(P2);                       // Clear phenotype data we allocated
          updateR(0);
        }
        difcormatrix[m][p] = chiSQ(ngenotypes, cors, nsamples);
        free(cors); free(nsamples);       // Clear correlation and samples data we allocated
      }
    }
//    if(verbose) info("Done marker %d/%d\n", m, genotypes.nmarkers);
    for(g = 0; g < ngenotypes; g++){      // Clear splits and pheno1
      free(splits[g].data);
    }
    freematrix((void*)pheno1, ngenotypes);
    free(splits);
  }
  return difcormatrix;
}

