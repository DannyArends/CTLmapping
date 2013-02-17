/******************************************************************//**
 * \file Rctl/src/sort.h
 * \brief Definition of helper functions related to sorting of 2D vectors
 *
 * <i>Copyright (c) 2010-2013</i>GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#ifdef __cplusplus
  extern "C" {
#endif
  #ifndef __SORT_H__
    #define __SORT_H__

    #include "ctl.h"

    int c_cmp(const void *a, const void *b);
    int d_cmp(const void *a, const void *b);
    int i_cmp(const void *a, const void *b);

  #endif //__SORT_H__
#ifdef __cplusplus
  }
#endif
