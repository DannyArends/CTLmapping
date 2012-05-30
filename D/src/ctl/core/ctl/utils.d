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

size_t[] which(int[] marker,int type = 0){
  size_t[] indices;
  indices.reserve(marker.length);
  for(size_t i=0; i < marker.length; i++){ if(marker[i] == type) indices ~= i; }
  return indices;
}

double[] get(double[] phenotype, size_t[] r){
  double[] ph;
  ph.reserve(r.length);
  foreach(size_t e;r){ ph ~= phenotype[e]; }
  return ph;
}

double mysign(double v){ if(v >=0) return 1.0; return -1.0; }
