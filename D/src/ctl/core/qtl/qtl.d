/******************************************************************//**
 * \file ctl/core/qtl/qtl.d
 * \brief Single marker QTL mapping (when qtlHD isn't found)
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.qtl.qtl;

import std.stdio, std.math, std.datetime, std.conv;
import ctl.core.array.matrix, ctl.core.array.ranges;
import ctl.core.qtl.utils, ctl.core.qtl.regression;
import ctl.core.analysis, ctl.io.terminal;

class SingleQTL : Analysis{
  override double[][] analyse(in int[][] genotypes,in double[][] phenotypes, in int[] geno_cov = [], bool verbose = true){
    if(verbose) MSG("Starting QTL mapping");
    SysTime stime = Clock.currTime();
    double[][] lodmatrix = newmatrix!double(phenotypes.length, genotypes.length, 0);
    if(verbose) write(" ");
    size_t npheno = phenotypes.length;
    for(size_t p=0; p < npheno; p++){
      for(size_t m=0; m < genotypes.length; m++){
        double[] w = newvector!double(phenotypes[0].length, 1.0);
        int[] nm = newvector!int(1, 1);
        lodmatrix[p][m] = multipleregression(createdesignmatrix(genotypes, m, geno_cov), phenotypes[p], w, nm, false);
      }
      if((p % max!uint((npheno/20),to!uint(1))) == 0 && verbose) write(".");
      stdout.flush();
    }
    if(verbose) writeln();
    if(verbose) MSG("QTL mapping done (%s secs)",(Clock.currTime()-stime).total!"seconds"());
    return lodmatrix;
  }
}
double[][] createdesignmatrix(in int[][] genotypes, int marker, in int[] geno_cov = [], bool intercept = true){
  double[][] dm;
  dm.length = genotypes[0].length;
  size_t ncols = 1 + geno_cov.length + cast(size_t)intercept;
  for(size_t v=0; v < ncols; v++){
    for(size_t i=0; i < genotypes[0].length; i++){
      dm[i].length = 1 + geno_cov.length + cast(size_t)intercept;
      if(intercept && v==0){
        dm[i][v] = 1.0;
      }else{
        if(v==(dm[i].length-1)){
          dm[i][v] = cast(double) genotypes[marker][i];
        }else{
          uint cov = v - cast(uint) intercept;
          dm[i][v] = cast(double) genotypes[geno_cov[cov]][i];
        }
      }
    }
  }
  return dm;
}
