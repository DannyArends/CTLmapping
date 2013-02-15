#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __VECTOR_H__
    #define __VECTOR_H__

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    #include <float.h>
    #include "ctl.h"

    //Custom length index vector
    typedef struct{
      int* data;
      size_t  nelements;
    } clvector;

    char*    newcvector(size_t dim);
    double*  newdvector(size_t dim);
    int*     newivector(size_t dim);

    char*    addtocvector(const char*   v, size_t dim, char n);    
    double*  addtodvector(const double* v, size_t dim, double n);
    int*     addtoivector(const int*    v, size_t dim, int n);
    
    void     printcvector(const char*   v, size_t dim);
    void     printdvector(const double* v, size_t dim);
    void     printivector(const int*    v, size_t dim);

    int*     randomizeivector(int* idx, size_t dim);

    /* Which integer elements in v are equal to type, returns a clvector of indices */
    inline clvector which(const int* v, size_t dim, int type){
      size_t i  = 0;
      clvector clv;
      for(i = 0; i < dim; i++){ 
        if(v[i] == type){
          clv.data = addtoivector(clv.data, clv.nelements, i);
          clv.nelements++;
        }
      }
      return clv;
    }
    
    /* Get the double elements in v, specified by the indexes in the clvector idxs */
    inline double* get(const double* v, clvector idxs){
      size_t i;
      double* v1 = newdvector(idxs.nelements);
      for(i = 0; i < idxs.nelements; i++){
        v1[i] = v[idxs.data[i]]; 
      }
      return v1;
    }

    inline int in(const clvector vector, int val){
      size_t i;
      for(i =0; i< vector.nelements; i++){
        if(val == vector.data[i] && val != -999) return 1;
      }
      return 0;
    }
  #endif //__VECTOR_H__
#ifdef __cplusplus
  }
#endif
