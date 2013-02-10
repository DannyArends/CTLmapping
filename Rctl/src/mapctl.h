#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __MAPCTL_H__
    #define __MAPCTL_H__

    #include "ctl.h"
    #include "correlation.h"
    #include "permutation.h"    
    #include "sort.h"
    void     updateR(int flush);
    void     R_mapctl(int* nind, int* nmar, int* nphe, int* geno, double* pheno, int* p, int *nperms, int* a, int* g, int* permt, double* dcor, double* perms, double* res, int* verb);
    double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, int alpha, int gamma, int nperms);
    double** diffcor(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, int alpha, int gamma);
		void call_from_d(int x);

  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

