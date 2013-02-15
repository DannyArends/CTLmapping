#include "structs.h"
#include "ctlio.h"
#include "mapctl.h"
#include <getopt.h>

/* Print the help, listing all possible command line switches */
void printhelp(){
  printf("Usage:\n");
  printf(" mapctl -g<genotype_file> -p<phenotype_file>\n\n");
  printf(" -p<FILE>   Input file with phenotype data (Default: phenotypes.csv)\n");
  printf(" -g<FILE>   Input file with genotype data (Default: genotypes.csv)\n");
  printf(" -n<N>      # of permutations (Default: 100)\n");
  printf(" -h         Shows this help\n");
}

/* Main function of the command line tool */
int main(int argc, char **argv){
  printf("Correlated Trait Locus (CTL) mapping\n");
  printf("(c) 2012 GBIC, written by Danny Arends\n");
  printf("Number of command line arguments passed: %d\n", (argc-1));
  char*  genofilename = "../D/test/data/genotypes.csv";
  char*  phenofilename= "../D/test/data/phenotypes.csv";
  size_t nperms = 100;
  size_t alpha  = 1;
  size_t beta   = 1;
  char   ch;
  srand(time(NULL));
  while((ch = getopt(argc, argv, "g:p:n:a:b:h")) != -1){
    switch(ch){
      case 'g': genofilename  = optarg;  break;
      case 'p': phenofilename = optarg; break;
      case 'n': nperms        = atoi(optarg);  break;
      case 'a': alpha         = atoi(optarg);  break;
      case 'b': beta          = atoi(optarg);  break;
      case 'h': printhelp(); return 0;  break;
      default: break;
    }
  }
  
  char* phenocontent    = getfilecontent(phenofilename);
  char* genocontent     = getfilecontent(genofilename);
  Phenotypes phenotypes = tophenotypes(phenocontent);
  Genotypes genotypes   = togenotypes(genocontent);
  free(phenocontent);
  free(genocontent);
  if(phenotypes.nindividuals != genotypes.nindividuals){
    printf("Individuals doesn't match between genotypes and phenotypes");
    return -1;
  }else{
    clvector genoenc = getGenotypes(genotypes);
    info("Num genotypes: %d\n", genoenc.nelements);
    size_t p = 0;
    for(p = 0; p < phenotypes.nphenotypes;p++){
      double** ctls = mapctl(phenotypes, genotypes, p, genoenc.nelements, genoenc.data, alpha, beta, nperms);
      writeout(ctls, p, genotypes.nmarkers, phenotypes.nphenotypes);
      freematrix((void**)ctls, genotypes.nmarkers);
    }
    freematrix((void**)phenotypes.data, phenotypes.nphenotypes);
    freematrix((void**)genotypes.data, genotypes.nmarkers);
  }
  printf("All done. Thank you for using mapctl\n");
  printf("Please cite: CTL mapping - Journal - Arends et al. [2013]\n");
  return 0;
}

