/**********************************************************************
 * src/ctl/core/stats/pxpmatrix.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written Feb, 2012
 **********************************************************************/ 
module ctl.core.stats.pxpmatrix;

import std.stdio;
import std.conv;
import std.array;
import std.datetime;

double[][] topxpmatrix(double[][][] ctlmmatrix, bool verbose = true){
  SysTime stime = Clock.currTime();
  double[][] m;

  if(verbose) writeln(" - PxP matrix: ",(Clock.currTime()-stime).total!"msecs"()," msecs");  
  return m;
}
