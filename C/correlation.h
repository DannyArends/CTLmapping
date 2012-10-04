#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CORRELATION_H__
    #define __CORRELATION_H__
    
    #include "structs.h"
    
    double correlation(double* x, double* y, size_t dim);
  
  #endif //__CORRELATION_H__
#ifdef __cplusplus
  }
#endif