/******************************************************************//**
 * \file ctl/mapctl.d
 * \brief Function holding the main mapctl workflow
 *
 * <i>Copyright (c) 2012</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written May, 2011<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
import std.stdio, std.conv, std.getopt, std.file, std.datetime, core.memory, std.string;
import std.algorithm, std.array;
import ctl.core.array.matrix, ctl.core.ctl.utils;
import ctl.io.reader, ctl.io.csv.write, ctl.io.csv.parse;

extern (C){
  struct Phenotypes{
    double** data;
    size_t   nphenotypes;
    size_t   nindividuals;
  };

  struct clvector{
    int*   data;
    size_t nelements;
  }

  struct Genotypes{
    int**    data;
    size_t   nmarkers;
    size_t   nindividuals;
  };

  clvector* getGenotypes(const Genotypes geno, bool verbose);

  double**  ctleffects(const Phenotypes phenotypes, const Genotypes genotypes, size_t phenotype, 
                       clvector* genoenc, int alpha, int beta, bool verbose);

  double*   permute(const Phenotypes phe, const Genotypes geno, size_t p, clvector* genoenc, 
                    int a, int b, size_t np, bool verbose);

  double**  toLODexact(double** scores, clvector* genoenc, size_t nmar, size_t nphe);
  double**  toLOD(double** scores, double* permutations, size_t nmar, size_t nphe, size_t nperms);
}

void main(string[] args){
  SysTime stime = Clock.currTime();
  writeln("Correlated Trait Locus (CTL) mapping in D");
  writeln("(c) 2012 written by Danny Arends in the D programming language");
	bool help      = false;
  bool verbose   = false;
  bool overwrite = false;
	uint alpha     = 1;
	uint beta      = 1;
  bool doperms   = false;
	uint nperms    = 100;
  string outdir  = "ctlout";
  string genofilename  = "./test/data/genotypes.csv";
  string phenofilename = "./test/data/phenotypes.csv";
	string format        = "csv";

  getopt(args, "help|h"     , &help
             , "verbose|v"  , &verbose
             , "redo|r"     , &overwrite
             , "alpha"      , &alpha
             , "beta"       , &beta
             , "nperms|n"   , &nperms
             , "out|o"      , &outdir
             , "geno|g"     , &genofilename
             , "pheno|p"    , &phenofilename
             , "format|f"   , &format);

  Reader ireader  = initialize(format);

  if(overwrite) writefln("Overwriting files in on");
  if(verbose){ writefln("Verbose mode on"); writefln("Output saved to: " ~ outdir); }
  if(!help){
    writefln("Start loading input files (%s, %s)",  phenofilename, genofilename);
    double[][]  phenotypes = ireader.loadphenotypes(phenofilename);
    string[]    phenonames = ireader.loadphenonames(phenofilename);
    int[][]     genotypes  = ireader.loadgenotypes(genofilename);
    if(verbose) writefln("Dataset: %s geno- and %s phenotypes", genotypes.length, phenotypes.length);
    if(genotypes.length == 0) abort("No genotypes loaded, analysis aborted"); 
    if(phenotypes.length == 0) abort("No phenotypes loaded, analysis aborted");
    if(verbose) writefln("Dataset: %s and %s individuals", genotypes[0].length, phenotypes[0].length);
    if(genotypes[0].length != phenotypes[0].length){
      abort(xformat("Mismatch between individuals %s != %s", genotypes[0].length, phenotypes[0].length));
    }
    Genotypes  geno;
    Phenotypes pheno;

    geno.data = toPP!int(genotypes);                             // Fill the genotype structure for C
    geno.nmarkers     = genotypes.length;
    geno.nindividuals = genotypes[0].length;

    pheno.data = toPP!double(phenotypes);                        // Fill the phenotype structure for C
    pheno.nphenotypes  = phenotypes.length;
    pheno.nindividuals = phenotypes[0].length;

    if(!exists(outdir)) mkdir(outdir);
    clvector* genoenc = getGenotypes(geno, false);
  
    for(size_t p = 0; p < phenotypes.length; p++){               // Main CTL mapping loop
      if(verbose) writef("- Phenotype %s: Mapping",p);
      string fnctl    = outdir ~ "/ctls"~to!string(p)~".out";
      string fnlods   = outdir ~ "/lods"~to!string(p)~".out";

      double**   ctlp = ctleffects(pheno, geno, p, genoenc, alpha, beta, false);
      double[][] ctls = translate(fromPP(ctlp, pheno.nphenotypes, geno.nmarkers));
      writeFile(ctls, fnctl, null, overwrite, verbose);
      double[][] lods;
      if(doperms){
        write(", Permutation");
        double*   permp = permute(pheno, geno, p, genoenc, alpha, beta, nperms, false);

        writeln(", toLOD");
        double**   lodp = toLOD(ctlp, permp, geno.nmarkers, pheno.nphenotypes, nperms);
        lods = translate(fromPP(lodp, pheno.nphenotypes, geno.nmarkers));
        
      }else{
        double**   lodp = toLODexact(ctlp, genoenc, geno.nmarkers, pheno.nphenotypes);
        lods = translate(fromPP(lodp, pheno.nphenotypes, geno.nmarkers));
        writeFile(lods, fnlods, null, overwrite, verbose);
      }
    }
    writefln("\nDone after %s seconds. Thank you for using mapctl", (Clock.currTime()-stime).total!"seconds"());
    writeln("Please cite: CTL mapping - Journal - Arends et al. [2013]");
  }
}

