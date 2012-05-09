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
  if(output[($-1)]=='/') output = output[0..($-1)];
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
    //We need an output path default to ./ ??
    if(!exists(output)) mkdirRecurse(output);

    //Start by mapping all QTL
    if(needanalysis(output ~ "/qtls.txt",overwrite)){
      Analysis a = getanalysis(settings);
      double[][] result  = a.analyse(genotypes, phenotypes, [], verbose);
      writeFile(result, output ~ "/qtls.txt", overwrite, verbose);
    }else{ writeln("\n[SKIP] QTL mapping"); }
    //Do the CTL
    double[][][] ctlmmatrix;
    for(uint p=0; p < phenotypes.length; p++){
      if(verbose) write("-Phenotype ",p);
      double[][] score;
      double[][] perms;
      
      if(needanalysis(output ~ "/ctl"~to!string(p)~".txt",overwrite)){
        score = mapping(phenotypes,  genotypes, p, verbose);
        writeFile(translate(score),  output ~ "/ctl"~to!string(p)~".txt", overwrite, verbose);
      }else{ writeln("\n[SKIP] CTL mapping"); }
      
      if(needanalysis(output ~ "/perms"~to!string(p)~".txt",overwrite)){
        perms = permutation(phenotypes, genotypes, p, settings.getInt("--nperms"), verbose);
        writeFile(translate(perms),  output ~ "/perms"~to!string(p)~".txt", overwrite, verbose);
      }else{ writeln("\n[SKIP] Permutations"); }
      
      if(needanalysis(output ~ "/lodscores"~to!string(p)~".txt",overwrite)){
        ctlmmatrix  ~= tolod(score, perms, verbose);
        writeFile(translate(ctlmmatrix[p]),  output ~ "/lodscores"~to!string(p)~".txt", overwrite, verbose);
      }else{ writeln("\n[SKIP] LOD transformation"); }
    }
    writeln("\nmapCTL finished analysis took: ",(Clock.currTime()-stime).total!"seconds"()," seconds");
    writefln("Continue by starting R and loading the results:\n library(ctl)\n ctls <- load.ctl(\"%s\", \"%s\", \"%s\")\n plot(ctls)",input_g, input_p,output);
  }
}
