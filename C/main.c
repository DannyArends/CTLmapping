#include "structs.h"
#include "ctlio.h"
#include "correlation.h"
#include <getopt.h>
 
int main(int argc, char **argv){
  printf("Correlated Trait Locus (CTL) mapping\n");
  printf("(c) 2012 written by Danny Arends\n");
  printf("Number of command line arguments passed: %d\n", argc);
  char* genofilename = "../D/test/data/genotypes.csv";
  char* phenofilename= "../D/test/data/phenotypes.csv";
  char   ch;
  size_t cnt;
  while((ch = getopt(argc, argv, "g:p:")) != -1){
    switch(ch){
      case 'g':
        genofilename = optarg;
      break;
      case 'p':
        phenofilename = optarg;
      break;
      default: break;
    }
  }
  
  char* phenocontent  = getfilecontent(phenofilename);
  double** phenotypes = todmatrix(phenocontent);
  char* genocontent   = getfilecontent(genofilename);
  int** genotypes     = toimatrix(genocontent);
  return 0;
}
