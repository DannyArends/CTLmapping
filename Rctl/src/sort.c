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
  return ceil(*ia  - *ib);
}

/* qsort int comparison function */
int i_cmp(const void *a, const void *b){
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;
  return *ia  - *ib;
}
