/**********************************************************************
 * src/ctl/core/qtl/nrc.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.core.qtl.nrc;

import std.stdio;
import std.math;
import std.mathspecial;

double inverseF(int df1, int df2, double alfa, bool verbose = false){
  double prob=0.0, minF=0.0, maxF=100.0, halfway=50.0, absdiff=1.0;
 int count=0;
  while ((absdiff>0.001)&&(count<100)){
    count++;
    halfway= (maxF+minF)/2.0;
    prob= betaIncomplete(df2/2.0,df1/2.0,df2/(df2+df1*halfway));
    if (prob<alfa) maxF= halfway;
    else minF= halfway;
    absdiff= fabs(prob-alfa);
  }
  if(verbose) writeln("prob=" , prob , "; alfa=" , alfa);
  return halfway;
}
