/******************************************************************//**
 * \file ctl/core/ctl/mapping.d
 * \brief Map CTL for a single phenotypes
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.ctl.mapping;

import std.stdio, std.math, std.datetime;
import ctl.core.array.matrix, ctl.io.terminal;
import ctl.core.ctl.utils, ctl.core.array.search;
import ctl.core.stats.correlation;

int[][] getEncodings(in int[][] genotypes, bool verbose = true){
  SysTime stime = Clock.currTime();
  int[][] r;
  for(size_t m=0; m<genotypes.length; m++){
    int[] d;
    foreach(int geno; genotypes[m]){ if(geno != int.max && !searchArray(d,geno)) d ~= geno; }
    r ~= d;
  }
  if(verbose) MSG("Scan of encoded genotype took: (%s msecs)\n",(Clock.currTime()-stime).total!"msecs"());
  return r;
}

double[][] mapping(in double[][] phenotypes, in int[][] genotypes, in int[][] encodings, uint phenotype = 1, bool verbose = true){
  assert(phenotype < phenotypes.length);
  SysTime stime = Clock.currTime();
  double[][] difcormatrix = newmatrix!double(genotypes.length, phenotypes.length, 0);
  for(size_t m=0; m<genotypes.length; m++){
    if(encodings[m].length >= 2){
      for(size_t x=0;x<(encodings[m].length-1);x++){
      for(size_t y=(x+1);y<encodings[m].length;y++){
      size_t[] ind_aa  = which(genotypes[m],encodings[m][x]);
      size_t[] ind_bb  = which(genotypes[m],encodings[m][y]);
      if(ind_aa.length < 2 || ind_bb.length < 2){
        for(size_t p=0; p<phenotypes.length; p++){ difcormatrix[m][p] = 0.0; }
      }else{
        double[] pheno_aa = get(phenotypes[phenotype],ind_aa);
        double[] pheno_bb = get(phenotypes[phenotype],ind_bb);
        for(size_t p=0; p<phenotypes.length; p++){
          double cor_aa = correlation!double(pheno_aa, get(phenotypes[p],ind_aa));
          double cor_bb = correlation!double(pheno_bb, get(phenotypes[p],ind_bb));
          if(isnan(cor_aa) || isnan(cor_bb)){
            difcormatrix[m][p] = 0.0;
          }else{
            difcormatrix[m][p] = pow(cor_aa - cor_bb, 2);
          }
        }
      }
      }
      }
    }
  }
  if(verbose) MSG("CTL mapping took: (%s msecs)",(Clock.currTime()-stime).total!"msecs"());
  return difcormatrix;
}
