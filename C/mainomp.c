#include <omp.h>
#include <cstdio>

int main(){
  const int nt = omp_get_max_threads();
  printf("openMP with %d threads", nt);

  #pragma omp parallel
  {
   printf("Hello world by %d", omp_get_thread_num());
  }
}

