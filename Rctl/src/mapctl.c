#include "mapctl.h"

/* Function to 'update' R, checks user input and can flushes console */
void updateR(int flush){
  #ifdef USING_R
    R_CheckUserInterrupt();
    if(flush) R_FlushConsole();
  #endif
}

/* R interface to perform a CTL scan and permutations on phenotype 'phenotype' */
void R_mapctl(int* nind, int* nmar, int* nphe, int* ngeno, int* geno, double* pheno, int* genoenc,
                      int* p, int *nperms, int* a, int* b, int* permt, double* dcor, 
                      double* perms, double* res, int* verb){
  int nindividuals  = (int)(*nind);
  int nmarkers      = (int)(*nmar);
  int nphenotypes   = (int)(*nphe);
  int ngenotypes    = (int)(*ngeno);
  int phenotype     = (int)(*p);
  int npermutations = (int)(*nperms);
  int alpha         = (int)(*a);
  int beta          = (int)(*b);
  int permtype      = (int)(*permt);
  int verbose       = (int)(*verb);

  Phenotypes phenotypes;
  Genotypes  genotypes;
  size_t     i,ph,perm;
  double**   dcors;
  double**   ctls;
  
  phenotypes.data = asdmatrix(nphenotypes, nindividuals, pheno);
  phenotypes.nphenotypes = nphenotypes;
  phenotypes.nindividuals = nindividuals;
  
  genotypes.data = asimatrix(nmarkers, nindividuals, geno);
  genotypes.nmarkers = nmarkers;
  genotypes.nindividuals = nindividuals;

  if(verbose) info("Phenotype %d: Mapping", (phenotype+1));  
  updateR(1);
  dcors = ctleffects(phenotypes, genotypes, phenotype, ngenotypes, genoenc, alpha, beta);
  for(i=0; i < (nphenotypes*nmarkers); i++){
    int m = i % nmarkers; int p = i / nmarkers;
    dcor[i] = dcors[m][p];
  }

  if(permtype){
    if(verbose) info(", RW permutation");
    updateR(1);
    double** permutations = permuteRW(phenotypes, genotypes, phenotype, ngenotypes, genoenc, alpha, beta, npermutations, 0);
    for(ph=0; ph < (nphenotypes); ph++){ // SEND PERMUTATIONS TO R
      for(perm=0; perm < (npermutations); perm++){
        perms[(ph*npermutations)+perm] = permutations[ph][perm];
      }
    }
    if(verbose) info(", toLOD\n");
    updateR(1);
    ctls = toLODRW(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
    freematrix((void**)permutations, nphenotypes);
  }else{
    if(verbose) info(", Permutation");
    updateR(1);
    double* permutations = permute(phenotypes, genotypes, phenotype, ngenotypes, genoenc, alpha, beta, npermutations, 0);
    for(i=0; i < npermutations; i++){ // SEND PERMUTATIONS TO R
      perms[i] = permutations[i];
    }
    if(verbose) info(", toLOD\n");
    updateR(1);
    ctls = toLOD(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
    free(permutations);
  }
  for(i=0; i < (nphenotypes*nmarkers); i++){
    int m = i % nmarkers;
    int p = i / nmarkers;
    res[i] = ctls[m][p];
  }
  freematrix((void**)dcors, genotypes.nmarkers);
  freematrix((void**)ctls, genotypes.nmarkers);
  return;
}

/* Perform a CTL scan and permutations on phenotype 'phenotype' */
double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, size_t ngenotypes, int* genoenc, int alpha, int beta, int nperms){
  info("Phenotype %d: Mapping", (phenotype+1));
  double** scores = ctleffects(phenotypes, genotypes, phenotype, ngenotypes, genoenc, alpha, beta);
  if((alpha == 1 && beta == 1)){
    info("\n");
    return scores;
  }
  info(", Permutation");
  double* permutations = permute(phenotypes, genotypes, phenotype, ngenotypes, genoenc, alpha, beta, nperms, 0);
  info(", toLOD\n");
  double** ctls = toLOD(scores, permutations, genotypes.nmarkers, phenotypes.nphenotypes, nperms);

  freematrix((void**)scores, genotypes.nmarkers);
  free(permutations);
  return ctls;
}

/* Calculate the difference in correlation matrix for phenotype 'phenotype' */
double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, size_t ngenotypes, int* genoenc, int alpha, int beta){
  size_t g, m, p,debug = 0;
  double** difcormatrix = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);
  for(p = 0; p < phenotypes.nphenotypes; p++){  
    if(p!=phenotype){
    for(m = 0; m < genotypes.nmarkers; m++){
      if(debug) info("Search and allocation marker %d\n", p);
      double* cors  = newdvector(ngenotypes);
      int* nsamples = newivector(ngenotypes);
      for(g = 0; g < ngenotypes; g++){
        clvector ind = which(genotypes.data[m], phenotypes.nindividuals, genoenc[g]);
        double* P1   = get(phenotypes.data[phenotype], ind);
        double* P2   = get(phenotypes.data[p], ind);
        cors[g]      = correlation(P1, P2, ind.nelements);
        nsamples[g]  = ind.nelements;
        free(P1); free(P2);       // Clear the data we allocated
        free(ind.data);
        updateR(0);
      }
      difcormatrix[m][p] = chiSQ(ngenotypes, cors, nsamples);
      free(cors); free(nsamples);       // Clear the data we allocated
      if(debug) info("Marker done\n");
    }
    }
  }
  return difcormatrix;
}


