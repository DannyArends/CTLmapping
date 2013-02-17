/******************************************************************//**
 * \file Rctl/src/rmapctl.c
 * \brief Implementation of the interfaces to R, and the updateR helper function
 *
 * <i>Copyright (c) 2010-2013</i>GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "rmapctl.h"

/* Function to 'update' R, checks user input and can flushes console */
void updateR(bool flush){
  #ifdef USING_R
    R_CheckUserInterrupt();
    if(flush) R_FlushConsole();
  #endif
}

/* R interface to perform a CTL scan and permutations on phenotype 'phenotype' */
void R_mapctl(int* nind, int* nmar, int* nphe, int* geno, double* pheno, int* p, 
              int *nperms, int* a, int* b, int* permt, double* dcor,  double* perms, double* res, int* verb){

  int nindividuals  = (int)(*nind);
  int nmarkers      = (int)(*nmar);
  int nphenotypes   = (int)(*nphe);
//  int ngenotypes    = (int)(*ngeno);
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

  clvector* genoenc = getGenotypes(genotypes, false);

  if(verbose) info("Phenotype %d: Mapping", (phenotype+1));  
  updateR(1);
  dcors = ctleffects(phenotypes, genotypes, phenotype, genoenc, alpha, beta, verbose);

  for(i=0; i < (nphenotypes*nmarkers); i++){    // Send scores to R
    int m = i % nmarkers; int p = i / nmarkers;
    dcor[i] = dcors[m][p];
  }

  if(permtype == 1){
    if(verbose) info(", Full permutation");
    updateR(1);
    double* permutations = permute(phenotypes, genotypes, phenotype, genoenc, 
                                   alpha, beta, npermutations, false);
    for(i=0; i < npermutations; i++){           // Send permutations to R
      perms[i] = permutations[i];
    }
    if(verbose) info(", toLOD\n");
    updateR(1);
    ctls = toLOD(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
    free(permutations);
  }else if(permtype == 2){
    if(verbose) info(", Pairwise permutation");
    updateR(1);
    double** permutations = permuteRW(phenotypes, genotypes, phenotype, genoenc, alpha, beta, npermutations, false);
    for(ph=0; ph < (nphenotypes); ph++){         // Send permutations to R
      for(perm=0; perm < (npermutations); perm++){
        perms[(ph*npermutations)+perm] = permutations[ph][perm];
      }
    }
    if(verbose) info(", toLOD\n");
    ctls = toLODRW(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
    freematrix((void**)permutations, nphenotypes);
  }else{
    if(verbose) info(", toLOD\n");
    updateR(1);
    ctls = toLODexact(dcors, genoenc, genotypes.nmarkers, phenotypes.nphenotypes);
  }
  for(i=0; i < (nphenotypes*nmarkers); i++){
    int m = i % nmarkers;
    int p = i / nmarkers;
    res[i] = ctls[m][p];
  }
  for(i = 0; i < nmarkers; i++){ free(genoenc[i].data); }
  free(genoenc);
  freematrix((void**)dcors, genotypes.nmarkers);
  freematrix((void**)ctls, genotypes.nmarkers);
  return;
}

