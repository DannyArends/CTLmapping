/**********************************************************************
 * src/qcl/core/qcl/singleqcl.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module qcl.core.qcl.singleqcl;

import std.stdio;
import std.math;
import std.datetime;

import qcl.core.qcl.utils;
import qcl.core.stats.correlation;

double[][] singleqcl(double[][] phenotypes, int[][] genotypes, uint phenotype = 1, bool verbose = true){
  assert(phenotype < phenotypes.length);
  SysTime stime = Clock.currTime();
  double[][] difcormatrix;
  for(uint m=0; m<genotypes.length; m++){
    double[] sqcordifference;
    uint[] ind_aa  = which(genotypes[m],0);
    uint[] ind_bb  = which(genotypes[m],1);
    double[] pheno_aa = get(phenotypes[phenotype],ind_aa);
    double[] pheno_bb = get(phenotypes[phenotype],ind_bb);
    for(uint p=0; p<phenotypes.length; p++){
      double cor_aa =correlation!double(pheno_aa, get(phenotypes[p],ind_aa));
      double cor_bb =correlation!double(pheno_bb, get(phenotypes[p],ind_bb));
      sqcordifference ~= ((mysign(cor_aa)*pow(cor_aa,2)) - (mysign(cor_bb)*pow(cor_bb,2)));
    }
    difcormatrix ~= sqcordifference;
  }
  if(verbose) writeln(" - Differential correlations: ",(Clock.currTime()-stime).total!"msecs"()," msecs");
  return difcormatrix;
}
