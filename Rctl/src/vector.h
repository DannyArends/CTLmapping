/******************************************************************//**
 * \file Rctl/src/vector.h
 * \brief Definition of functions related to 2D vectors
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __VECTOR_H__
    #define __VECTOR_H__

    #include "ctl.h"

    /** Representation of an custom length integer vector. */
    typedef struct{
      int*    data;               /*!< Integer data*  */
      size_t  nelements;          /*!< Number of elements in data  */
    } clvector;

    /** Allocate a new character vector of length dim */
    char*    newcvector(size_t dim);
    /** Allocate a new double vector of length dim */
    double*  newdvector(size_t dim);
    /** Allocate a new integer vector of length dim */
    int*     newivector(size_t dim);
    /** Allocate a new clvector of length dim */
    clvector newclvector(size_t dim);

    /** Add a new character element n to vector v */
    char*    addtocvector(char*   v, size_t dim, char n);
    /** Ads a new double element n to vector v */
    double*  addtodvector(double* v, size_t dim, double n);
    /** Add a new integer element n to vector v */
    int*     addtoivector(int*    v, size_t dim, int n);
    
    /** Print dim elements from a character vector to the output */
    void     printcvector(const char*   v, size_t dim);
    /** Print dim elements from a double vector to the output */
    void     printdvector(const double* v, size_t dim);
    /** Print dim elements from a integer vector to the output */
    void     printivector(const int*    v, size_t dim);
    /** Print a custom length vector to the output */
    void printclvector(const clvector v);

    /** Get the indices from v where v[i]==e */
    clvector which(const int* v, size_t dim, int e);
    /** Get by index a new double* containing elements from v */
    double*  get(const double* v, clvector idxs);
    /** Return true when an element in v matches e. */
    bool     in(const clvector v, int e);
    /** Get the maximum value in vector v */
    double   vectormax(double* v, size_t dim);
    /** Fisher-Yates algorithm to generate a random-range */
    int*     randomizeivector(int* idx, size_t dim);

    /** Double* to ranked int*. This is a possibly slow implementation 
     *  to get the ranks in a double vector while we propagates our 
     *  missing value (-999).<br>
     *  e.g.: 5      7    8   3   5
     *  to:   1.5    2    3   1   1.5 */
    double*  torank(double* v, size_t dim);

  #endif //__VECTOR_H__
#ifdef __cplusplus
  }
#endif

