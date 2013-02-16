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
  size_t nperms  = 100;
  size_t alpha   = 1;
  size_t beta    = 1;
  bool   doperms = false;
  bool   verbose = false;
  char   ch;

  srand(time(NULL));
  while((ch = getopt(argc, argv, "g:p:n:a:b:hdv")) != -1){
    switch(ch){
      case 'g': genofilename  = optarg;  break;
      case 'p': phenofilename = optarg; break;
      case 'd': doperms       = true;  break;
      case 'v': verbose       = true;  break;
      case 'n': nperms        = atoi(optarg);  break;
      case 'h': printhelp(); return 0;  break;
      default: break;
    }
  }
  
  char* phenocontent    = getfilecontent(phenofilename);
  char* genocontent     = getfilecontent(genofilename);
  Phenotypes phenotypes = tophenotypes(phenocontent);
  free(phenocontent);
  Genotypes genotypes   = togenotypes(genocontent);
  free(genocontent);

  if(phenotypes.nindividuals != genotypes.nindividuals){
    printf("Individuals doesn't match between genotypes and phenotypes");
    return -1;
  }else{
    clvector genoenc = getGenotypes(genotypes);
    info("Num genotypes: %d\n", genoenc.nelements);
    size_t phenotype  = 0;
    size_t ngenos = genoenc.nelements;
    for(phenotype = 0; phenotype < phenotypes.nphenotypes; phenotype++){
      info("Phenotype %d: Mapping", (phenotype+1));

      double** ctls;
      double*  perms;
      double** scores = ctleffects(phenotypes, genotypes, phenotype, ngenos, genoenc.data, alpha, beta, verbose);
      if(!doperms){
        info(", toLOD\n");  // Exact calculation can be used
        ctls = toLODexact(scores, ngenos, genotypes.nmarkers, phenotypes.nphenotypes);
      }else{
        info(", Permutation");
        fflush(stdout);
        perms = permute(phenotypes, genotypes, phenotype, ngenos, genoenc.data, alpha, beta, nperms, false);
        info(", toLOD\n");
        ctls = toLOD(scores, perms, genotypes.nmarkers, phenotypes.nphenotypes, nperms);
        free(perms);
      }
      writeout(ctls, phenotype, genotypes.nmarkers, phenotypes.nphenotypes);
      writesummary(ctls, scores, phenotype, genotypes.nmarkers, phenotypes.nphenotypes, 0.05);
      freematrix((void**)scores, genotypes.nmarkers);
      freematrix((void**)ctls, genotypes.nmarkers);
    }
    info("%.8f\n",1-chiSQtoP(1,80.798));
    freematrix((void**)phenotypes.data, phenotypes.nphenotypes);
    freematrix((void**)genotypes.data, genotypes.nmarkers);
  }
  printf("All done. Thank you for using mapctl\n");
  printf("Please cite: CTL mapping - Journal - Arends et al. [2013]\n");
  return 0;
}

