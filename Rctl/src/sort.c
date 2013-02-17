/******************************************************************//**
 * \file Rctl/src/sort.c
 * \brief Implementation of helper functions related to sorting of 2D vectors
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "sort.h"

/* qsort char / char* comparison function */
int c_cmp(const void *a, const void *b){
  const char **ia = (const char **)a;
  const char **ib = (const char **)b;
  return strcmp(*ia, *ib);
}

/* qsort double comparison function */
int d_cmp(const void *a, const void *b){
  const double *ia = (const double *)a;
  const double *ib = (const double *)b;
  if(*ia < *ib) return -1;
  else if(*ia > *ib) return 1;
  return 0;
}

/* qsort int comparison function */
int i_cmp(const void *a, const void *b){
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;
  return *ia  - *ib;
}
