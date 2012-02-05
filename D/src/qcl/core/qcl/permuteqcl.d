/**********************************************************************
 * src/qcl/core/qcl/singleqcl.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module qcl.core.qcl.permuteqcl;

import std.stdio;
import std.math;
import std.datetime;
import std.random;

import qcl.core.array.ranges;
import qcl.core.array.matrix;
import qcl.core.qcl.singleqcl;
import qcl.core.stats.basic;

T[][] permute(T)(T[][] genotypes){
  T[][] newgeno = newmatrix!T(genotypes.length,genotypes[0].length);
  uint cnt=0;
  foreach(e; randomCover(dorange(0, (genotypes[0].length-1)),Random(unpredictableSeed))){
    for(uint g=0;g<genotypes.length;g++){
      newgeno[g][e] = genotypes[g][cnt];
    }
    cnt++;
  }
  return newgeno;
}

double[][] permuteqcl(double[][] phenotypes, int[][] genotypes, uint phenotype = 1, uint permutations = 100, bool verbose = true){
  assert(phenotype < phenotypes.length);
  SysTime stime = Clock.currTime();
  double[][] permutationmatrix;
  if(verbose) write(" ");
  for(uint p=0; p<permutations; p++){
    double[][] pqclm = singleqcl(phenotypes, permute!int(genotypes), phenotype, false);
    permutationmatrix ~= doMatrixMax!double(pqclm);
    if(p % 3 == 0 && verbose){write("."); stdout.flush();}
  }
  if(verbose) writeln("\n - Permutations: ",(Clock.currTime()-stime).total!"seconds"(),"secs");
  return permutationmatrix;
}
