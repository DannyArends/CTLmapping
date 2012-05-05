/**********************************************************************
 * src/ctl/mapctl.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified May, 2012
 * first written May, 2011
 **********************************************************************/ 
import std.stdio;
import std.math;
import std.conv;
import std.datetime;

import ctl.core.array.matrix;
import ctl.core.stats.basic;
import ctl.core.stats.tolod;
import ctl.core.stats.pxpmatrix;
import ctl.core.analysis;
import ctl.core.ctl.mapping;
import ctl.core.ctl.permutation;
import ctl.io.reader;
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
  Reader r = initialize(settings);
  double[][]  phenotypes = r.loadphenotypes(settings.getString("--phenotypes"));
  int[][]     genotypes  = r.loadgenotypes(settings.getString("--genotypes"));
  if(!settings.displayHelp()){
    bool verbose = settings.getBool("--verbose");
    bool overwrite = settings.getBool("--overwrite");

    if(overwrite) writeln("[overwrite] Overwriting files in on");
    if(verbose){
      writefln("[verbose] on");
      writefln("[info] %s geno- and %s phenotypes", genotypes.length, phenotypes.length);
      writefln("[info] Measurements on %s and %s individuals\n", genotypes[0].length, phenotypes[0].length);
    }
    assert(genotypes[0].length == phenotypes[0].length, "mismatch between individuals");
    //Start by mapping all QTL
    Analysis a = getanalysis(settings);
    double[][] result  = a.analyse(genotypes, phenotypes, [], verbose);
    writeFile(result,  "test/output/qtls.txt", overwrite, verbose);
    //Do the CTL
    double[][][] ctlmmatrix;
    for(uint p=0; p < phenotypes.length; p++){
      if(verbose) write("-Phenotype ",p);
      double[][] score = mapping(phenotypes,  genotypes, p, verbose);
      double[][] perms = permutation(phenotypes, genotypes, p, settings.getInt("--nperms"), verbose);
      ctlmmatrix  ~= tolod(score, perms, verbose);
      writeFile(translate(ctlmmatrix[p]),  "test/output/lodscores"~to!string(p)~".txt", overwrite, verbose);
    }
    double[][] pxpmatrix = topxpmatrix(ctlmmatrix);
    writeFile(pxpmatrix,  "test/output/pxpmatrix.txt", overwrite, verbose);    
    writeln("\nmapCTL finished analysis took: ",(Clock.currTime()-stime).total!"seconds"()," seconds");
  }
}

