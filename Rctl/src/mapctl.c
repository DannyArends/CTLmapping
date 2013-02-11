#include "mapctl.h"

/* Function to 'update' R, checks user input and can flushes console */
void updateR(int flush){
  #ifdef USING_R
    R_CheckUserInterrupt();
    if(flush) R_FlushConsole();
  #endif
}

/* R interface to perform a CTL scan and permutations on phenotype 'phenotype' */
void R_mapctl(int* nind, int* nmar, int* nphe, int* geno, double* pheno, int* p, int *nperms, int* a, int* b, int* permt, double* dcor, double* perms, double* res, int* verb){
  int nindividuals  = (int)(*nind);
  int nmarkers      = (int)(*nmar);
  int nphenotypes   = (int)(*nphe);
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
  dcors = ctleffects(phenotypes, genotypes, phenotype, alpha, beta);
  for(i=0; i < (nphenotypes*nmarkers); i++){
    int m = i % nmarkers; int p = i / nmarkers;
    dcor[i] = dcors[m][p];
  }

  if(!(alpha == 1 && beta == 1)){
    if(permtype){
      if(verbose) info(", RW permutation");
      updateR(1);
      double** permutations = permuteRW(phenotypes, genotypes, phenotype, alpha, beta, npermutations, 0);
      for(ph=0; ph < (nphenotypes); ph++){ // SEND PERMUTATIONS TO R
        for(perm=0; perm < (npermutations); perm++){
          perms[(ph*npermutations)+perm] = permutations[ph][perm];
        }
      }
      if(verbose) info(", toLOD\n");
      updateR(1);
      ctls = toLODRW(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
    }else{
      if(verbose) info(", Permutation");
      updateR(1);
      double* permutations = permute(phenotypes, genotypes, phenotype, alpha, beta, npermutations, 0);
      for(i=0; i < npermutations; i++){ // SEND PERMUTATIONS TO R
        perms[i] = permutations[i];
      }
      if(verbose) info(", toLOD\n");
      updateR(1);
      ctls = toLOD(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
    }
  }else{ // alpha == 1 && beta == 1
    ctls = dcors;
    info("\n");
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
double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, int alpha, int beta, int nperms){
  info("Phenotype %d:", (phenotype+1));
  double** scores = ctleffects(phenotypes, genotypes, phenotype, alpha, beta);
  if((alpha == 1 && beta == 1)){
    info("\n");
    return scores;
  }
  info(", Permutation");
  double* permutations = permute(phenotypes, genotypes, phenotype, alpha, beta, nperms, 0);
  info(", toLOD\n");
  double** ctls = toLOD(scores, permutations, genotypes.nmarkers, phenotypes.nphenotypes, nperms);

  freematrix((void**)scores, genotypes.nmarkers);
  free(permutations);
  return ctls;
}

/* Calculate the standard error */
double stderror(size_t df1, size_t df2){
  return(sqrt((1.0 / ((double)(df1-3)) + (1.0 / (double)(df2-3)))));
}

/* Transform a correlation coeficient into a Zscore */
double zscore(double cor){ return(.5*log((1.0 + cor)/(1.0 - cor))); }


double ctleff(double* phe1, double* phe2, int* m, int nind, int alpha, int beta, int doZ){
  clvector indAA  = which(m, nind, 0);
  clvector indBB  = which(m, nind, 1);
  double* phe1AA = get(phe1, indAA);
  double* phe1BB = get(phe1, indBB);
  double* phe2AA = get(phe2, indAA);
  double* phe2BB = get(phe2, indBB);
  double  cAA = correlation(phe1AA, phe2AA, indAA.nelements);
  double  cBB = correlation(phe1BB, phe2BB, indBB.nelements);
  if(doZ && alpha == 1 && beta == 1){
    return (zscore(cAA) - zscore(cBB)) / stderror(indAA.nelements, indBB.nelements);
  }
  return pow(0.5*(pow(cAA, alpha) - pow(cBB, alpha)), beta);
}

/* Calculate the difference in correlation matrix for phenotype 'phenotype' */
double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, int alpha, int beta){
  size_t m,p,debug = 0;
  double** difcormatrix = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);
  for(m = 0; m < genotypes.nmarkers; m++){
    if(debug) info("Search and allocation marker %d\n", p);
    clvector ind_aa  = which(genotypes.data[m], phenotypes.nindividuals, 0);
    clvector ind_bb  = which(genotypes.data[m], phenotypes.nindividuals, 1);
    double* pheno_aa1 = get(phenotypes.data[phenotype],ind_aa);
    double* pheno_bb1 = get(phenotypes.data[phenotype],ind_bb);
    for(p = 0; p < phenotypes.nphenotypes; p++){
      if(p!=phenotype){
        if(debug) info("Search and allocation phenotype %d", p);
        double* pheno_aa2 = get(phenotypes.data[p],ind_aa);
        double* pheno_bb2 = get(phenotypes.data[p],ind_bb);
        if(debug) info(", done");
        double cor_aa = correlation(pheno_aa1, pheno_aa2, ind_aa.nelements);
        double cor_bb = correlation(pheno_bb1, pheno_bb2, ind_bb.nelements);
        if(debug) info(", correlation done");
        if(alpha == 1 && beta == 1){// DEFAULT No permutations, just exact calculations
          difcormatrix[m][p] = (zscore(cor_aa) - zscore(cor_bb)) / stderror(ind_aa.nelements, ind_bb.nelements);
        }else{
          difcormatrix[m][p] = pow(0.5*(pow(cor_aa, alpha) - pow(cor_bb, alpha)), beta);
        }
        if(debug) info(", cleanup\n");
        free(pheno_aa2);
        free(pheno_bb2);
        updateR(0);
      }
    }
    if(debug) info("Marker done\n");
    free(pheno_aa1);
    free(pheno_bb1);
    free(ind_aa.data);
    free(ind_bb.data);
  }
  return difcormatrix;
}

