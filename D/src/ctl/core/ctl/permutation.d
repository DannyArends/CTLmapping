/**********************************************************************
 * src/ctl/core/ctl/permute.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.core.ctl.permutation;

import std.stdio;
import std.math;
import std.datetime;
import std.random;

import ctl.core.array.ranges;
import ctl.core.array.matrix;
import ctl.core.ctl.mapping;
import ctl.core.stats.basic;

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

double[][] permutation(double[][] phenotypes, int[][] genotypes, uint phenotype = 1, uint permutations = 100, bool verbose = true){
  assert(phenotype < phenotypes.length);
  SysTime stime = Clock.currTime();
  double[][] permutationmatrix;
  if(verbose) write(" ");
  for(uint p=0; p<permutations; p++){
    double[][] perm_m = mapping(phenotypes, permute!int(genotypes), phenotype, false);
    permutationmatrix ~= doMatrixMax!double(perm_m);
    if(p % 3 == 0 && verbose){write("."); stdout.flush();}
  }
  if(verbose) writeln("\n - Permutations: ",(Clock.currTime()-stime).total!"seconds"(),"secs");
  return permutationmatrix;
}
