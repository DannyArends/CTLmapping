/******************************************************************//**
 * \file ctl/core/stats/tolod.d
 * \brief Transformation of a CTL matrix into a LOD score matrix
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.stats.tolod;
 
import std.stdio, std.math, std.conv, std.array, std.datetime;
import ctl.core.array.ranges, ctl.core.array.matrix, ctl.core.array.search;
import ctl.io.terminal, ctl.core.stats.basic, ctl.core.ctl.utils;

double[][] tolod(double[][] scores, double[] permlist, size_t[] significant, string[] phenonames, bool verbose = true){
  SysTime stime = Clock.currTime();
  double[][] m = [];
  foreach(size_t s; significant){
    double[] profile;
    for(size_t c=0;c<scores[0].length;c++){
      profile ~= getEstimate(permlist,abs(scores[s][c]));
    }
    m ~= profile;
  }
  MSG("Above threshold: %s traits ( %s )", m.length, get(phenonames,significant));
  if(verbose) MSG("LOD transformation: (%s msecs)",(Clock.currTime()-stime).total!"msecs"());  
  return m;
}

size_t[] getSignificant(double[][] scores, double critval){
  size_t[] indices = null;
  for(size_t r=0;r<scores.length;r++){
    if(max!double(scores[r]) >= critval){
      indices ~= r;
    }
  }
  return indices;
}

double getEstimate(in double[] permlist, double score){
  return abs(log10(1 - (to!double(getIndex(permlist, score))/to!double(permlist.length+1))));
}

double getCutoff(in double[] permlist, double minlod){
  double fdr = pow(10,-minlod);
  double critval = permlist[to!int((1.0-fdr)*permlist.length)];
  double estimate = getEstimate(permlist,critval);
  MSG("FDR %s|%s, %s, %s|%s",pow(10,-estimate), fdr, critval, getEstimate(permlist,critval), minlod);
  return critval;
}
