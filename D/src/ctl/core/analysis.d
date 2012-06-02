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

import std.stdio, std.file, std.math, std.datetime,std.conv;
import ctl.core.qtl.qtl;
import ctl.core.stats.basic;
import ctl.io.cmdline.parse, ctl.io.terminal;
import ctl.core.array.matrix, ctl.core.array.ranges;
import ctl.core.ctl.utils;

bool needanalysis(string filename, bool overwrite = false){
  if(!exists(filename)) return true;
  if(overwrite) return true;
  return false;
}

Analysis getanalysis(string analysis, CTLsettings settings){
  switch(analysis){
    case "qtl":
      if(settings.getBool("--qtl")) return new SingleQTL();
    break;
    case "effect":
      if(settings.getBool("--effect")) return new EffectScan();
    break;
    default:
    break;
  }
  return new NullAnalysis();
}

abstract class Analysis{
  double[][] analyse(in int[][] genotypes, in double[][] phenotypes, in int[] geno_cov, bool verbose);
}

class NullAnalysis : Analysis{
  override double[][] analyse(in int[][] genotypes, in double[][] phenotypes, in int[] geno_cov = [], bool verbose = true){
    double[][] m = null;
    return m;
  }
}

class EffectScan : Analysis{
  override double[][] analyse(in int[][] genotypes, in double[][] phenotypes, in int[] geno_cov = [], bool verbose = true){
    double[][] effects = newmatrix!double(genotypes.length, phenotypes.length, 0);
    if(verbose) MSG("Starting effect scan");
    SysTime stime = Clock.currTime();
    size_t ngeno = genotypes.length;
    if(verbose) write(" ");
    for(size_t m=0; m<ngeno; m++){
      size_t[] ind_aa  = which(genotypes[m],0);
      size_t[] ind_bb  = which(genotypes[m],1);
      for(size_t p=0; p < phenotypes.length; p++){
        double m_aa = doMean!double(get(phenotypes[p],ind_aa));
        double m_bb = doMean!double(get(phenotypes[p],ind_bb));
        effects[m][p] = m_aa - m_bb;
      }
      if((m % max!uint((ngeno/20),to!uint(1))) == 0 && verbose) write(".");
      stdout.flush();
    }
    writeln();
    if(verbose) MSG("Effect scan took: (%s msecs)",(Clock.currTime()-stime).total!"msecs"());
    return effects;
  }
}
