/******************************************************************//**
 * \file C/main.c
 * \brief Main entry point for the C standalone version
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "structs.h"
#include "vector.h"
#include "ctlio.h"
#include "mapctl.h"
#include "permutation.h"
#include <getopt.h>

/** Print the help, listing all possible command line switches. */
void printhelp(){
  printf("Usage:\n");
  printf(" mapctl -g<genotype_file> -p<phenotype_file>\n\n");
  printf(" -p<FILE>   Input file with phenotype data (Default: phenotypes.csv)\n");
  printf(" -g<FILE>   Input file with genotype data (Default: genotypes.csv)\n");
  printf(" -o<FILE>   Name of the output file (Default: summary.txt)\n");
  printf(" -t<N>      Significance threshold (0..1) (Default: 0.01)\n");
  printf(" -f         Save all the results to file \n");
  printf(" -d         When specified permutation are performed \n");
  printf(" -n<N>      # of permutations (Default: 100)\n");
  printf(" -h         Shows this help\n");
}

/** Main function of the command line tool. */
int main(int argc, char **argv){
  printf("Correlated Trait Locus (CTL) mapping\n");
  printf("(c) 2012-2020 GBIC - RUG & HU-Berlin, written by Danny Arends\n");
  printf("Number of command line arguments passed: %d\n", (argc-1));
  char *r_argv[] = { "R", "--silent" };
  Rf_initEmbeddedR(2, r_argv);
#ifdef TEST
  char*  genofilename  = "../D/test/data/genotypes.csv";
  char*  phenofilename = "../D/test/data/phenotypes.csv";
#else
  char*  genofilename  = "genotypes.csv";
  char*  phenofilename = "phenotypes.csv";
#endif
  char*  outfilename   = "summary.txt";
  size_t nperms        = 100;
  size_t alpha         = 1;
  size_t beta          = 1;
  bool   doperms       = false;
  bool   verbose       = false;
  bool   fulloutput    = false;
  bool   nooutput      = false;
  double threshold     = 0.01;
  char   ch;

  srand(time(NULL));
  while((ch = getopt(argc, argv, "g:p:n:a:b:t:hdvfs")) != -1){
    switch(ch){
      case 'g': genofilename  = optarg;  break;
      case 'p': phenofilename = optarg; break;
      case 'd': doperms       = true;  break;
      case 'f': fulloutput    = true;  break;
      case 's': nooutput      = true;  break;
      case 'v': verbose       = true;  break;
      case 't': threshold     = atof(optarg);  break;
      case 'n': nperms        = atoi(optarg);  break;
      case 'h': printhelp(); return 0;  break;
      default: break;
    }
  }
  
  char* phenocontent    = getFilecontent(phenofilename);
  char* genocontent     = getFilecontent(genofilename);
  Phenotypes phenotypes = toPhenotypes(phenocontent);
  free(phenocontent);
  Genotypes genotypes   = toGenotypes(genocontent);
  free(genocontent);

  if(phenotypes.nindividuals != genotypes.nindividuals){
    printf("Individuals doesn't match between genotypes and phenotypes");
    return -1;
  }else{


    size_t nmar = genotypes.nmarkers;
    size_t nphe = phenotypes.nphenotypes;
    size_t m;

    int rthreads = omp_get_max_threads();                         /* Requested number of threads */
    int nitems = nphe;                              /* Number of items todo */
    int npthread = ceil((float)nitems / (float)rthreads);         /* Number we do in every thread */


    clvector* genoenc = getGenotypes(genotypes, verbose);

    #pragma omp parallel
    {
    size_t phenotype  = 0;
    int tid = omp_get_thread_num();                                                 /* Obtain thread number */
    int start = npthread * tid;
    int stop = (int)fmin((float)(start + (int)npthread), (float)nitems);
    info("I am thread = %d/%d -> [%d,%d]\n", tid, rthreads, start, stop);
    for(size_t phenotype = start; phenotype < stop; phenotype++){
//      info("Phenotype %d: Mapping", (phenotype+1));

      double** ctls;
      double*  perms;
      double** scores = ctleffects(phenotypes, genotypes, phenotype, genoenc, 2, verbose);
      if(!doperms){
        info(", toLOD\n", "");  // Exact calculation can be used
        ctls = toLODexact(scores, genoenc, nmar, nphe, false);
      }else{
        info(", Permutation", "");
        perms = permute(phenotypes, genotypes, phenotype, genoenc, nperms, 1, false);
        info(", toLOD\n", "");
        ctls = toLOD(scores, perms, nmar, nphe, nperms);
        free(perms);
      }
      if(!nooutput){
        if(fulloutput) writeout(ctls, phenotype, nmar, nphe);
        writesummary(phenotypes, genotypes, outfilename, ctls, phenotype, nmar, nphe, genoenc, threshold);
      }  
      freematrix((void**)scores, genotypes.nmarkers);
      freematrix((void**)ctls, genotypes.nmarkers);
      //info("Phenotype %d done\n", (phenotype+1));
    }
    //info("Thread %d done\n", tid);
    } //OMP

    for(m = 0; m < nmar; m++){ free(genoenc[m].data); }
    free(genoenc);
    freematrix((void**)phenotypes.data, nphe);
    freematrix((void**)genotypes.data, nmar);
  }
  info("All done. Thank you for using mapctl\n", "");
  info("Please cite: CTL mapping - Journal - Arends et al. [2013]\n", "");
  Rf_endEmbeddedR(0);
  return 0;
}

