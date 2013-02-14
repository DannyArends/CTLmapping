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
    void     R_mapctl(int* nind, int* nmar, int* nphe, int* ngeno, int* geno, double* pheno, int* genoenc,
                      int* p, int *nperms, int* a, int* b, int* permt, double* dcor, 
                      double* perms, double* res, int* verb);

    double   ctleff(double* phe1, double* phe2, int* m, int nind, int alpha, int beta, int doZ);
    double   chiSQ(size_t nr, double* r, int* nsamples);

    double** mapctl(Phenotypes phenotypes, Genotypes genotypes, size_t phenotype, size_t ngenotypes, int* genoenc, int alpha, int beta, int nperms);
    double** ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, size_t ngenotypes, int* genoenc, int alpha, int beta);

  #endif //__MAPCTL_H__
#ifdef __cplusplus
  }
#endif

