/**********************************************************************
 * src/ctl/core/ctl/utils.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.core.ctl.utils;

import std.stdio;

import ctl.core.array.matrix;
import ctl.core.stats.correlation;

uint[] which(int[] marker,int type = 0){
  uint[] indices;
  indices.reserve(marker.length);
  for(uint i=0; i < marker.length; i++){
    if(marker[i] == type) indices ~= i; 
  }
  return indices;
}

double[] get(double[] phenotype, uint[] r){
  double[] ph;
  ph.reserve(r.length);
  foreach(uint e;r){
    ph ~= phenotype[e];
  }
  return ph;
}

double mysign(double v){
  if(v >=0) return 1.0;
  return -1.0;
}
