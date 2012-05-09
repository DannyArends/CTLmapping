/**********************************************************************
 * src/ctl/core/analysis.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
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
  double[][] analyse(int[][] genotypes, double[][] phenotypes, int[] geno_cov = [], bool verbose = true){
    double[][] m = newmatrix!double(1,1);
    return m;
  }
}
