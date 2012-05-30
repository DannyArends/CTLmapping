/******************************************************************//**
 * \file ctl/core/ctl/permutation.d
 * \brief CTL significance thresholds by permutation
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.ctl.permutation;

import std.stdio, std.math, std.datetime, std.random;
import ctl.core.array.ranges, ctl.io.terminal;
import ctl.core.array.matrix, ctl.core.ctl.mapping;
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
  permutationmatrix.reserve(permutations);
  if(verbose) write(" ");
  for(uint p=0; p<permutations; p++){
    double[][] perm_m = mapping(phenotypes, permute!int(genotypes), phenotype, false);
    permutationmatrix ~= doMatrixMax!double(perm_m);
    if(p % 3 == 0 && verbose){write("."); stdout.flush();}
  }
  if(verbose) writeln();
  if(verbose) MSG("Permutations took: (%s secs)", (Clock.currTime()-stime).total!"seconds"());
  return permutationmatrix;
}

