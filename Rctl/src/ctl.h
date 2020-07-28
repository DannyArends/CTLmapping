/******************************************************************//**
 * \file Rctl/src/ctl.h
 * \brief Global definitions and includes
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __CTL_H__
    #define __CTL_H__

    #include <stdbool.h>
    #include <string.h>

    #ifdef STANDALONE
      #include <stdio.h>
      #include <stdlib.h>
      #include <Rmath.h>
      #include <float.h>
      #define MISSING NAN
      #define CHECKNA isnan

      #define info(format, ...) { \
        printf(format, __VA_ARGS__); \
        fflush(stdout); }
      #define err(format, ...) { \
        printf(format, __VA_ARGS__); \
        exit(-1); }
    #else
      #define USING_R
      #include <R.h>
      #include <Rmath.h>
      #define MISSING R_NaN
      #define CHECKNA ISNA

      #define info(format, ...) { \
        Rprintf(format, __VA_ARGS__);}
      #define err(format, ...) { \
        error(format, __VA_ARGS__);}
      #endif
    
  #endif //__CTL_H__
#ifdef __cplusplus
  }
#endif

