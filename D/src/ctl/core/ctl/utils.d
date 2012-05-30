/******************************************************************//**
 * \file ctl/core/ctl/utils.d
 * \brief Utility classes for CTL mapping
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.core.ctl.utils;

import std.stdio;
import ctl.core.array.matrix;
import ctl.core.stats.correlation;

uint[] which(int[] marker,int type = 0){
  uint[] indices;
  indices.reserve(marker.length);
  for(uint i=0; i < marker.length; i++){ if(marker[i] == type) indices ~= i; }
  return indices;
}

double[] get(double[] phenotype, uint[] r){
  double[] ph;
  ph.reserve(r.length);
  foreach(uint e;r){ ph ~= phenotype[e]; }
  return ph;
}

double mysign(double v){ if(v >=0) return 1.0; return -1.0; }
