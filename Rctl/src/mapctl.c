#include "mapctl.h"

void updateR(int flush){
  #ifdef USING_R
    R_CheckUserInterrupt();
    if(flush) R_FlushConsole();
  #endif
}

void R_mapctl(int* nind, int* nmar, int* nphe, int* geno, double* pheno, int* p, int *nperms, int* a, int* g, int* permt, double* dcor, double* perms, double* res, int* verb){
  int nindividuals  = (int)(*nind);
  int nmarkers      = (int)(*nmar);
  int nphenotypes   = (int)(*nphe);
  int phenotype     = (int)(*p);
  int npermutations = (int)(*nperms);
  int alpha         = (int)(*a);
  int gamma         = (int)(*g);
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
  dcors = diffcor(phenotypes, genotypes, phenotype, alpha, gamma);
  for(i=0; i < (nphenotypes*nmarkers); i++){
    int m = i % nmarkers; int p = i / nmarkers;
    dcor[i] = floorf(dcors[m][p] * 1000 + 0.5) / 1000;
  }

  if(permtype){
    if(verbose) info(", Row-wise permutation");
    updateR(1);
    double** permutations = permuteRowWise(phenotypes, genotypes, phenotype, alpha, gamma, npermutations, 0);
    for(ph=0; ph < (nphenotypes); ph++){ // SEND PERMUTATIONS TO R
      for(perm=0; perm < (npermutations); perm++){
        perms[(ph*npermutations)+perm] = permutations[ph][perm];
      }
    }
    if(verbose) info(", toLOD\n");
    updateR(1);
    ctls = toLODRowWise(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
  }else{
    if(verbose) info(", Permutation");
    updateR(1);
    double* permutations = permute(phenotypes, genotypes, phenotype, alpha, gamma, npermutations, 0);
    for(i=0; i < npermutations; i++){ // SEND PERMUTATIONS TO R
      perms[i] = permutations[i];
    }
    if(verbose) info(", toLOD\n");
    updateR(1);
    ctls = toLOD(dcors, permutations, genotypes.nmarkers, phenotypes.nphenotypes, npermutations);
  }
  for(i=0; i < (nphenotypes*nmarkers); i++){
    int m = i % nmarkers;
    int p = i / nmarkers;
    res[i] = floorf(ctls[m][p] * 1000 + 0.5) / 1000;
  }
  freematrix((void**)dcors, genotypes.nmarkers);
  freematrix((void**)ctls, genotypes.nmarkers);
  return;
}

double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, int alpha, int gamma, int nperms){
  info("Phenotype %d: Mapping", (phenotype+1));
  double** dcorscores = diffcor(phenotypes, genotypes, phenotype, alpha, gamma);
  info(", Permutation");
  double* permutations = permute(phenotypes, genotypes, phenotype, alpha, gamma, nperms, 0);
  info(", toLOD\n");
  double** ctls = toLOD(dcorscores, permutations, genotypes.nmarkers, phenotypes.nphenotypes, nperms);
  freematrix((void**)dcorscores, genotypes.nmarkers);
  free(permutations);
  return ctls;
}

double** diffcor(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, int alpha, int gamma){
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
        if(alpha == 1 && gamma == 1){// DEFAULT No permutations, just exact calculations
          double z1 = .5*log((1.0 + cor_aa)/(1.0 - cor_aa));
          double z2 = .5*log((1.0 + cor_bb)/(1.0 - cor_bb));
          double se = sqrt((1.0 / (double)(ind_aa.nelements-3)) + (1.0 / (double)(ind_bb.nelements-3)));
          difcormatrix[m][p] = (z1 - z2) / se;
          //info("Score: %f %f %f -> %f\n",z1,z2,se,difcormatrix[m][p]);
        }else{
          difcormatrix[m][p] = pow(.5*(pow(cor_aa, alpha) - pow(cor_bb, alpha)), gamma);
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
