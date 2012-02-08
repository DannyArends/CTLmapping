/**********************************************************************
 * src/ctl/mapctl.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/ 
import std.stdio;
import std.math;
import std.conv;
import std.datetime;

import ctl.core.array.matrix;
import ctl.core.stats.basic;
import ctl.core.stats.tolod;
import ctl.core.ctl.mapping;
import ctl.core.ctl.permutation;
import ctl.io.csv.parse;
import ctl.io.cmdline.parse;
import ctl.io.csv.write;

void printHelp(){
  writeln("Documentation: http://www.dannyarends.nl/CTL/index.html");
  writeln("Usage:");
  writeln("mapCTL [permutations]\n");
  writeln(" [permutations] Number of permutations");
}

void main(string[] args){
  SysTime stime = Clock.currTime();
  writeln("mapCTL: Correlated Trait Locus (CTL) mapping in D");
  writeln("(c) 2012 written by Danny Arends in the D programming language");
  CTLsettings settings   = parseCmd(args);
  double[][]  phenotypes = parseFile!double(settings.getString("--phenotypes"));
  int[][]     genotypes  = parseFile!int(settings.getString("--genotypes"));
  bool        verbose    = settings.getBool("--verbose");
  if(verbose) writefln("%s geno- and %s phenotypes on (%s/%s) individuals\n", genotypes.length, phenotypes.length, genotypes[0].length, phenotypes[0].length);
  assert(genotypes[0].length == phenotypes[0].length);
  for(uint p=0; p < phenotypes.length; p++){
    if(verbose) writeln("-Phenotype ",p);
    double[][] score = mapping(phenotypes,  genotypes, p, verbose);
    double[][] perms = permutation(phenotypes, genotypes, p, settings.getInt("--nperms"), verbose);
    double[][] lod   = tolod(score, perms, verbose);
    writeFile(translate(lod),  "test/output/lodscores"~to!string(p)~".txt");
  }
  writeln("\nmapCTL finished analysis took: ",(Clock.currTime()-stime).total!"seconds"()," seconds");
}
