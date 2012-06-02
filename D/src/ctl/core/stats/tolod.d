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
import ctl.io.terminal, ctl.core.array.ranges, ctl.core.array.matrix, ctl.core.array.search;

double[][] tolod(double[][] scores, double[][] permutations, bool verbose = true){
  SysTime stime = Clock.currTime();
  double[] permlist = unlist(absmatrix(permutations));
  double[][] m = newmatrix!double(scores.length,scores[0].length);
  permlist.sort;
  for(size_t r=0;r<scores.length;r++){
    for(size_t c=0;c<scores[0].length;c++){
      auto index = getIndex(permlist,abs(scores[r][c]));
      double estimate = 1.0;
      if(index < permlist.length){
        estimate -= (to!double(index)/to!double(permlist.length));
      }else{ /*TODO GDP on top 10% of scores */      
        estimate -= (to!double(permlist.length-1)/to!double(permlist.length));
      }
      m[r][c] = abs(log10(estimate));
    }
  }
  if(verbose) MSG("LOD transformation: (%s msecs)",(Clock.currTime()-stime).total!"msecs"());  
  return m;
}


double[][] tolod2(double[][] scores, double[][] permutations, double minlod, bool verbose = true){
  SysTime stime = Clock.currTime();
  double[] permlist = unlist(absmatrix(permutations));
  double[][] m = [];
  permlist.sort;
  double percentage = pow(10,-minlod);
  MSG("perc -> %s",percentage);
  double critval = permlist[to!int((1.0-percentage)*permlist.length)];
  MSG("%s",critval);
  for(size_t r=0;r<scores.length;r++){
    if(max!double(scores[r]) >= critval){
      double[] profile;
      for(size_t c=0;c<scores.length;c++){
        size_t index = getIndex(permlist,abs(scores[r][c]));
        double estimate = 1.0;
        if(index < permlist.length){
          estimate -= (to!double(index)/to!double(permlist.length));
        }else{ /*TODO GDP on top 10% of scores */      
          estimate -= (to!double(permlist.length-1)/to!double(permlist.length));
        }
        profile ~= abs(log10(estimate));
      }
      m ~= profile;
    }
  }
  if(verbose) MSG("LOD transformation: (%s msecs)",(Clock.currTime()-stime).total!"msecs"());  
  return m;
}


