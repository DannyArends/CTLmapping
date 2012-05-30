/**********************************************************************
 * src/ctl/core/stats/tolod.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/ 
module ctl.core.stats.tolod;
 
import std.stdio, std.math, std.conv, std.array, std.datetime;
import ctl.io.terminal, ctl.core.array.matrix, ctl.core.array.search;

double[][] tolod(double[][] scores, double[][] permutations, bool verbose = true){
  SysTime stime = Clock.currTime();
  double[] permlist = unlist(absmatrix(permutations));
  double[][] m = newmatrix!double(scores.length,scores[0].length);
  permlist.sort;
  for(uint r=0;r<scores.length;r++){
    for(uint c=0;c<scores[0].length;c++){
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

