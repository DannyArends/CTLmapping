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

import std.stdio, std.math, std.datetime, std.random, std.conv;
import ctl.core.array.ranges, ctl.io.terminal, core.memory;
import ctl.core.array.matrix, ctl.core.ctl.mapping;
import ctl.core.stats.basic;

T[][] permute(T)(in T[][] genotypes){
  T[][] newgeno = newmatrix!T(genotypes.length,genotypes[0].length);
  uint cnt=0;
  foreach(e; randomCover(dorange(0, genotypes[0].length),Random(unpredictableSeed))){
    for(size_t g = 0; g < genotypes.length; g++){
      newgeno[g][e] = genotypes[g][cnt];
    }
    cnt++;
  }
  return newgeno;
}

double[][] permutation(in double[][] phenotypes, in int[][] genotypes, size_t phenotype = 1, uint permutations = 100, bool verbose = true){
  assert(phenotype < phenotypes.length);
  SysTime stime = Clock.currTime();
  double[][] perms;
  double[][] perm_t;
  if(verbose) write(" ");
  for(size_t p=0; p < permutations; p++){
    perm_t = mapping(phenotypes, permute!int(genotypes), phenotype, false);
    perms ~= doMatrixMax!double(perm_t);
    if((p % max!uint((permutations/20),to!uint(1))) == 0 && verbose){write("."); stdout.flush();}
    GC.collect();
    GC.minimize();
  }
  if(verbose) writeln();
  if(verbose) MSG("Permutations took: (%s secs)", (Clock.currTime()-stime).total!"seconds"());
  return perms;
}
