#include "structs.h"
#include "ctlio.h"
#include "mapctl.h"
#include <getopt.h>
 
int main(int argc, char **argv){
  printf("Correlated Trait Locus (CTL) mapping\n");
  printf("(c) 2012 written by Danny Arends\n");
  printf("Number of command line arguments passed: %d\n", argc);
  char* genofilename = "../D/test/data/genotypes.csv";
  char* phenofilename= "../D/test/data/phenotypes.csv";
  char   ch;
  while((ch = getopt(argc, argv, "g:p:")) != -1){
    switch(ch){
      case 'g': genofilename = optarg;  break;
      case 'p': phenofilename = optarg; break;
      default: break;
    }
  }
  
  char* phenocontent    = getfilecontent(phenofilename);
  char* genocontent     = getfilecontent(genofilename);
  Phenotypes phenotypes = todmatrix(phenocontent);
  Genotypes genotypes   = toimatrix(genocontent);
  if(phenotypes.nindividuals != genotypes.nindividuals){
    printf("Individuals doesn't match between genotypes and phenotypes");
    return -1;
  }else{
    size_t p;
    for(p = 0; p < phenotypes.nphenotypes;p++){
      printf("Mapping phenotype %d\n",p);
      double** scores = ctlmapping(phenotypes, genotypes, p);
      //printdmatrix(scores,10,10);
    }
  }
  return 0;
}
