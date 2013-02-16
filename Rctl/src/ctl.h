#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTL_H__
    #define __CTL_H__

    #include <stdbool.h>

    #define MISSING -999

    #ifdef STANDALONE
      #define info(format, ...) { \
        printf(format, ## __VA_ARGS__);}
      #define err(format, ...) { \
        printf(format, ## __VA_ARGS__); \
        exit(-1);}
    #else
      #define USING_R
      #include <R.h>
      #define info(format, ...) { \
        Rprintf(format, ## __VA_ARGS__);}
      #define err(format, ...) { \
        error(format, ## __VA_ARGS__);}
      #endif
    
  #endif //__CTL_H__
#ifdef __cplusplus
  }
#endif

