#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __PERMUTATION_H__
    #define __PERMUTATION_H__

    #include "mapctl.h"
    #include "sort.h"
    
    double*   permute(const Phenotypes phe, const Genotypes geno, size_t p, int a, int b, size_t np, int verbose);
    double**  permuteRW(const Phenotypes phe, const Genotypes geno, size_t p, int a, int b, size_t np, int verbose);

    double**  toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms);
    double**  toLODRW(double** scores, double** permutations, size_t nmar, size_t nphe, size_t nperms);
  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

