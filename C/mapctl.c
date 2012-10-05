#include "mapctl.h"

double random(){
  return rand() / (RAND_MAX+1.0);
}

int* swap(int* idx, int i1, int i2){
  int t = idx[i2];
  idx[i2] = idx[i1];
  idx[i1] = t;
  return idx;
}

//Fisher-Yates random-range generation
int* randomizerange(int* idx, int max){
  if(max==0) return idx;
  return randomizerange(swap(idx, (int)(random()*(max-1)), max),(max-1));
}

Genotypes permutegenotypes(Genotypes genotypes){
  size_t m,i;
  size_t* idx = newivector(genotypes.nindividuals);
  for(i = 0; i < genotypes.nindividuals; i++){ idx[i] = i; }

  idx = randomizerange(idx, genotypes.nindividuals);

  int** newgenodata = newimatrix(genotypes.nmarkers, genotypes.nindividuals);

  for(m = 0; m < genotypes.nmarkers; m++){
    for(i = 0; i < genotypes.nindividuals; i++){
      newgenodata[m][i] = genotypes.data[m][idx[i]];
    }
  }
  freevector((void*)idx);
  Genotypes g = genotypes;
  g.data = newgenodata;
  return g;
}

double matrixmax(double** m, size_t rows, size_t cols){
  size_t r,c;
  double max = -DBL_MAX;
  for(r = 0; r < rows; r++){
    for(c = 0; c < cols; c++){
      if(m[r][c] > max) max = m[r][c];
    }
  }
  return max;
}

double* permutation(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, size_t nperms, int verbose){
  size_t p;
  double* scores = newdvector(nperms);
  for(p = 0; p < nperms; p++){
    Genotypes g = permutegenotypes(genotypes);
    double** ctls = ctlmapping(phenotypes, g, phenotype);
    scores[p] = matrixmax(ctls,genotypes.nmarkers,phenotypes.nphenotypes);
    freematrix((void**)ctls   , genotypes.nmarkers);
    freematrix((void**)g.data , genotypes.nmarkers);
    if(verbose) printf("Done with permutation %d\n",p);
  }
  return scores;
}

double** toLOD(double** scores, double* permutations){
  return scores;
}
 
double** ctlmapping(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype){
  size_t m,p;
  double** difcormatrix = newdmatrix(genotypes.nmarkers, phenotypes.nphenotypes);
  for(m = 0; m < genotypes.nmarkers; m++){
    IndexVector ind_aa  = which(genotypes.data[m], phenotypes.nindividuals, 0);
    IndexVector ind_bb  = which(genotypes.data[m], phenotypes.nindividuals, 1);
    double* pheno_aa1 = get(phenotypes.data[phenotype],ind_aa);
    double* pheno_bb1 = get(phenotypes.data[phenotype],ind_bb);
    for(p = 0; p < phenotypes.nphenotypes; p++){
      double* pheno_aa2 = get(phenotypes.data[p],ind_aa);
      double* pheno_bb2 = get(phenotypes.data[p],ind_bb);
      double cor_aa = correlation(pheno_aa1, pheno_aa2, phenotypes.nindividuals);
      double cor_bb = correlation(pheno_bb1, pheno_bb2, phenotypes.nindividuals);
      difcormatrix[m][p] = pow(cor_aa - cor_bb, 2);
      free(pheno_aa2);
      free(pheno_bb2);
    }
    free(pheno_aa1);
    free(pheno_bb1);
    free(ind_aa.data);
    free(ind_bb.data);
  }
  return difcormatrix;
}
