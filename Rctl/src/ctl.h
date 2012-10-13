#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTL_H__
    #define __CTL_H__

    #ifdef STANDALONE
      #pragma message("NOT USING R")
      #define info(format, ...) { \
        printf(format, ## __VA_ARGS__);}
    #else
      #define USING_R
      #pragma message("USING R")
      #include <R.h>
      #define info(format, ...) { \
        Rprintf(format, ## __VA_ARGS__);}
    #endif
    
  #endif //__CTL_H__
#ifdef __cplusplus
  }
#endif

