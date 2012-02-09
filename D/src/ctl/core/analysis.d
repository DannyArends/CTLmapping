/**********************************************************************
 * src/ctl/core/analysis.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.core.analysis;

version(QTL){
  import ctl.core.qtl.qtl;
}
import ctl.io.cmdline.parse;
import ctl.core.array.matrix;

Analysis getanalysis(CTLsettings settings){
  switch(settings.getString("--analysis")){
version(QTL){
    case "qtab":
     return new QTL();
    break;
}
    default:
    break;
  }
  return new NullAnalysis();
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
