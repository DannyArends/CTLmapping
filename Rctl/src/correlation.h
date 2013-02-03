#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CORRELATION_H__
    #define __CORRELATION_H__

    #include <math.h>
    #include "ctl.h"
    #include "structs.h"
    
    double correlation(double* x, double* y, size_t dim);
    double tstat(double cor, int dim);
  #endif //__CORRELATION_H__
#ifdef __cplusplus
  }
#endif
