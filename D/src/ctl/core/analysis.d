/******************************************************************//**
 * \file ctl/core/analysis.d
 * \brief Interface definition of a genetic analysis algorithm
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.analysis;

import std.file;
version(QTL){
  import ctl.core.qtl.qtl;
}
import ctl.io.cmdline.parse;
import ctl.core.array.matrix;

bool needanalysis(string filename, bool overwrite = false){
  if(!exists(filename)) return true;
  if(overwrite) return true;
  return false;
}

Analysis getanalysis(CTLsettings settings){
version(QTL){
  return new SingleQTL();
}else{
  return new NullAnalysis();
}
}

abstract class Analysis{
  double[][] analyse(int[][] genotypes, double[][] phenotypes, int[] geno_cov, bool verbose);
}

class NullAnalysis : Analysis{
  override double[][] analyse(int[][] genotypes, double[][] phenotypes, int[] geno_cov = [], bool verbose = true){
    double[][] m = newmatrix!double(1,1,0);
    return m;
  }
}
