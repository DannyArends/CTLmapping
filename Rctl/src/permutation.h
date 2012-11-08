#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __PERMUTATION_H__
    #define __PERMUTATION_H__

    #include "mapctl.h"
    #include "sort.h"
    
    double*   permute(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, int alpha, int gamma, size_t nperms, int verbose);
    double**  permuteRowWise(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, int alpha, int gamma, size_t nperms, int verbose);

    double**  toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms);
    double**  toLODRowWise(double** scores, double** permutations, size_t nmar, size_t nphe, size_t nperms);
  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

