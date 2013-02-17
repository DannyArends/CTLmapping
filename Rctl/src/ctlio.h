#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTLIO_H__
    #define __CTLIO_H__

    #include "ctl.h"
    #include "rmapctl.h"
    #include "correlation.h"
    #include "structs.h"

    void   writeout(double** ctls, size_t phenotype, size_t nmar, size_t nphe);
    void   writesummary(const Phenotypes phenotypes, const Genotypes genotypes, const char* fn, double** ctls, 
                  size_t phenotype, size_t nmar, size_t nphe, clvector* genoenc, double cutoff);

    double* getCorrelations(const Phenotypes phenotypes, const Genotypes genotypes, size_t phe1, 
                            clvector genoenc, size_t mar, size_t phe2, bool verbose);

    size_t filesize(char* name);
    char*  getfilecontent(char* name);

  #endif //__CTLIO_H__
#ifdef __cplusplus
  }
#endif

