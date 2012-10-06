#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __MAPCTL_H__
    #define __MAPCTL_H__
    
    #include "correlation.h"
    #include "sort.h"
    
    double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, int nperms);
    double** diffcor(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype);
    double*  permutation(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, size_t nperms, int verbose);
    double** toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms);
  
  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

