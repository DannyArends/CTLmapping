    #ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __STRUCTS_H__
    #define __STRUCTS_H__

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <float.h>
    #include "ctl.h"
    #include "vector.h"
    #include "matrix.h"

    typedef struct{
      double** data;
      size_t   nphenotypes;
      size_t   nindividuals;
    } Phenotypes;

    typedef struct{
      int**    data;
      size_t   nmarkers;
      size_t   nindividuals;
    } Genotypes;

    typedef struct{
      double sum;
      double cor;
    } KahanAccumulator;

    /* KahanAccumulator - Helper function to create an empty Accumulator: see http://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    static inline KahanAccumulator createAccumulator(){
      KahanAccumulator t;
      t.sum = 0.0;
      t.cor = 0.0;
      return t;
    }

    /* KahanAccumulator - Using the Kahan summation algorithm: see http://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    static inline KahanAccumulator KahanSum(KahanAccumulator accumulation, double value){
      KahanAccumulator result;
      double y = value - accumulation.cor;
      double t = accumulation.sum + y;
      result.cor = (t - accumulation.sum) - y;
      result.sum = t;
      return result;
    }

    clvector getGenotypes(const Genotypes g);
    Phenotypes tophenotypes(const char* content);
    Genotypes  togenotypes(const char* content);

  #endif //__STRUCTS_H__
#ifdef __cplusplus
  }
#endif
