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
import std.file;
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

void main(string[] args){
  SysTime stime = Clock.currTime();
  writeln("mapCTL: Correlated Trait Locus (CTL) mapping in D");
  writeln("(c) 2012 written by Danny Arends in the D programming language");
  CTLsettings settings   = parseCmd(args);
  Reader r = initialize(settings);
  bool verbose    = settings.getBool("--verbose");
  bool overwrite  = settings.getBool("--overwrite");
  string output   = settings.getString("--output");
  string input_p  = settings.getString("--phenotypes");
  string input_g  = settings.getString("--genotypes");

  if(overwrite) writeln("[overwrite] Overwriting files in on");
  if(verbose){
    writefln("[verbose] on");
    writefln("[output] Output saved to: " ~ output);
  }
  if(!settings.displayHelp()){
    writefln("[load] Start with loading input files (%s, %s)",  input_p, input_g);
    double[][]  phenotypes = r.loadphenotypes(input_p);
    int[][]     genotypes  = r.loadgenotypes(input_g);
    if(verbose){
      writefln("[info] %s geno- and %s phenotypes", genotypes.length, phenotypes.length);
      writefln("[info] Measurements on %s and %s individuals\n", genotypes[0].length, phenotypes[0].length);
    }

    assert(genotypes[0].length == phenotypes[0].length, "mismatch between individuals");
    //Start by mapping all QTL
    Analysis a = getanalysis(settings);
    double[][] result  = a.analyse(genotypes, phenotypes, [], verbose);
    //We need an output path default to ./ ??
    if(!exists(output)) mkdirRecurse(output);
    writeFile(result, output ~ "/qtls.txt", overwrite, verbose);
    //Do the CTL
    double[][][] ctlmmatrix;
    for(uint p=0; p < phenotypes.length; p++){
      if(verbose) write("-Phenotype ",p);
      double[][] score = mapping(phenotypes,  genotypes, p, verbose);
      double[][] perms = permutation(phenotypes, genotypes, p, settings.getInt("--nperms"), verbose);
      ctlmmatrix  ~= tolod(score, perms, verbose);
      writeFile(translate(ctlmmatrix[p]),  output ~ "/lodscores"~to!string(p)~".txt", overwrite, verbose);
    }
    double[][] pxpmatrix = topxpmatrix(ctlmmatrix);
    writeFile(pxpmatrix, output ~ "pxpmatrix.txt", overwrite, verbose);
    writeln("\nmapCTL finished analysis took: ",(Clock.currTime()-stime).total!"seconds"()," seconds");
  }
}

