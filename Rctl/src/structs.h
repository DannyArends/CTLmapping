/******************************************************************//**
 * \file Rctl/src/structs.h
 * \brief Genotype and Phenotype structures definitions
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __STRUCTS_H__
    #define __STRUCTS_H__

    #include "ctl.h"
    #include "vector.h"
    #include "matrix.h"

    typedef struct{
      double** data;          /*! Holds the phenotype data as a double** */
      size_t   nphenotypes;   /*! Length of the first dimension of data */
      size_t   nindividuals;  /*! Length of the second dimension of data */
    } Phenotypes;

    typedef struct{
      int**    data;          /*! Holds the genotype data as a int** */
      size_t   nmarkers;      /*! Length of the first dimension of data */
      size_t   nindividuals;  /*! Length of the second dimension of data */
    } Genotypes;

    typedef struct{
      double sum;
      double cor;
    } KahanAccumulator;

    /** Create an empty Accumulator. 
     *  See http://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    inline KahanAccumulator createAccumulator(){
      KahanAccumulator t;
      t.sum = 0.0;
      t.cor = 0.0;
      return t;
    }

    /** Using the Kahan summation algorithm.
     *  See http://en.wikipedia.org/wiki/Kahan_summation_algorithm */
    inline KahanAccumulator KahanSum(KahanAccumulator accumulation, double value){
      KahanAccumulator result;
      double y = value - accumulation.cor;
      double t = accumulation.sum + y;
      result.cor = (t - accumulation.sum) - y;
      result.sum = t;
      return result;
    }

    /** Get a clvector* with per marker the genotype encodings. 
     *  The number of elements in the clvector is equal to the number of markers, 
     *  each clvector then holds the genotype encoding for that marker */
    clvector* getGenotypes(const Genotypes geno, bool verbose);

    /** Create the Phenotypes object.
     *  This function creates the Phenotypes struct from a \0 terminated character array */
    Phenotypes toPhenotypes(char* content);

    /** Create the Genotypes object.
     *  This function creates the Genotypes struct from a \0 terminated character array */
    Genotypes  toGenotypes(char* content);

  #endif //__STRUCTS_H__
#ifdef __cplusplus
  }
#endif
