#include "structs.h"
#include "ctlio.h"
#include "mapctl.h"
#include <getopt.h>
 
int main(int argc, char **argv){
  printf("Correlated Trait Locus (CTL) mapping\n");
  printf("(c) 2012 written by Danny Arends\n");
  printf("Number of command line arguments passed: %d\n", argc);
  char*  genofilename = "../D/test/data/genotypes.csv";
  char*  phenofilename= "../D/test/data/phenotypes.csv";
  size_t nperms = 10;
  char   ch;
  while((ch = getopt(argc, argv, "g:p:n:")) != -1){
    switch(ch){
      case 'g': genofilename = optarg;  break;
      case 'p': phenofilename = optarg; break;
      case 'n': nperms = atoi(optarg); break;
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
    size_t p = 0;
    for(p = 0; p < phenotypes.nphenotypes;p++){
      printf("Phenotype %d: Mapping",p);
      double** ctlscores = ctlmapping(phenotypes, genotypes, p);
      printf(", Permutation");
      double* permutations = permutation(phenotypes, genotypes, p, nperms, 0);
      printf(", toLOD\n");
      freematrix((void**)ctlscores, genotypes.nmarkers);
      free(permutations);
    }
    freematrix((void**)phenotypes.data, phenotypes.nphenotypes);
    freematrix((void**)genotypes.data, genotypes.nmarkers);
  }
  return 0;
}
