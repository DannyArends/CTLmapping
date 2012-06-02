/******************************************************************//**
 * \file ctl/core/qtl/utils.d
 * \brief Utility functions for QTL mapping
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.qtl.utils;

import std.stdio, std.math;
import ctl.core.array.matrix, ctl.core.qtl.nrc, ctl.core.qtl.LUdecomp;

double Lnormal(double residual, double variance){
  return exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
}

struct FITTED{
  double logL = 0.0;
  double[] Fy;
}

FITTED calculateloglikelihood(uint nsamples, in double[] residual, in double[] w, double variance, bool verbose = true){
  FITTED f;
  f.Fy   = newvector!double(nsamples, 0.0);
  double[] indL = newvector!double(nsamples, 0.0);

  for (uint i=0; i<nsamples; i++){
    f.Fy[i]  = Lnormal(residual[i],variance);
    indL[i]  += w[i] * f.Fy[i];
    f.logL   += log(indL[i]);
  }
  return f;
}

struct STATS{
  double   variance = 0.0;
  double[] fit;
  double[] residual;
}

STATS calculatestatistics(uint nvariables, uint nsamples, in double[][] xt, in double[] xtwy, in double[] y, in double[] w, bool verbose = true){
  STATS s;
  s.fit      = newvector!double(nsamples, 0.0);
  s.residual = newvector!double(nsamples, 0.0);

  for(uint i=0; i<nsamples; i++){
    s.fit[i]      = 0.0;
    s.residual[i] = 0.0;
    for(uint j=0; j<nvariables; j++){
      s.fit[i]     += xt[j][i] * xtwy[j];
    }
    s.residual[i]   = y[i]-s.fit[i];
    s.variance     += w[i]*pow(s.residual[i],2.0);
  }
  s.variance /= nsamples;
  return s;
}

double[] calculateparameters(uint nvariables, uint nsamples, in double[][] xt, in double[] w, in double[] y, bool verbose = true){
  int d=0;
  double xtwj;
  double[][] XtWX = newmatrix!double(nvariables, nvariables, 0.0);
  double[]   XtWY = newvector!double(nvariables, 0.0);
  int[]      indx = newvector!int(nvariables, 0);
  if(verbose) writefln(" - Calculating XtWX and XtWY");
  for(uint i=0; i<nsamples; i++){
    for(uint j=0; j<nvariables; j++){
      xtwj     = xt[j][i] * w[i];
      XtWY[j] += xtwj    * y[i];
      for(uint jj=0; jj<=j; jj++){
        XtWX[j][jj] += xtwj * xt[jj][i];
      }
    }
  }
  if(verbose) writeln(" - XtWX: ",XtWX);
  LUdecompose(XtWX, nvariables, indx, &d);
  LUsolve(XtWX, nvariables, indx, XtWY);
  
  return XtWY;
}
