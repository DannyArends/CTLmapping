/******************************************************************//**
 * \file ctl/io/cmdline/parse.d
 * \brief Functions to parse the commandline
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written May, 2011<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.io.cmdline.parse;

import std.stdio, std.string, std.conv;

struct S(T){
  string name;
  string description;
  T      value;
}

struct Scontainer{
  S!bool[]   booleans;
  S!uint[]   integers;
  S!double[] doubles;
  S!string[] strings;
}

struct CTLsettings{
  void load(Scontainer opt){
    this.opt = opt;
  }
  
  bool displayHelp(){
    if(getBool("--help")){
      foreach(S!bool s; opt.booleans){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      foreach(S!uint s; opt.integers){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      foreach(S!double s; opt.doubles){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      foreach(S!string s; opt.strings){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      return true;
    }
    return false;
  }
  
  uint getInt(string name){
    foreach(int cnt, S!uint s; opt.integers){
      if(name==s.name){ return s.value; }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return -1;
  }

  double getDouble(string name){
    foreach(int cnt, S!double s; opt.doubles){
      if(name==s.name){ return s.value; }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return double.nan;
  }
  
  bool getBool(string name){
    foreach(int cnt, S!bool s; opt.booleans){
      if(name==s.name){ return s.value; }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return false;
  }
  
  string getString(string name){
    foreach(int cnt, S!string s; opt.strings){
      if(name==s.name){ return s.value; }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return "";
  }
  Scontainer opt;
}

import std.getopt;

CTLsettings parseCmd(string[] args){
  CTLsettings settings;
  Scontainer opt;

  bool help = false;
  bool effect = false;
  bool qtl = false;
  bool verbose = true;
  bool redo = false;
  bool sp = false;
  uint nperms = 100;
  double minlod = 1.5;
  string phenotype_filename = "test/data/phenotypes.csv";
  string genotype_filename = "test/data/genotypes.csv";
  string output = "./";
  string file_format = "csv";

  getopt(args, "help|h", &help
             , "verbose|v", &verbose
             , "effect|e", &effect
             , "qtl|q", &qtl
             , "sp|s", &sp
             , "output|o", &output
             , "redo|r", &redo
             , "nperms|n", &nperms
             , "minlod", &minlod
             , "phenotypes|p", &phenotype_filename
             , "genotypes|g", &genotype_filename
             , "format|f", &file_format);

  opt.booleans ~= S!bool("--help"         ,"Show the help file", help);
  opt.booleans ~= S!bool("--verbose"      ,"Verbose mode", verbose);
  opt.booleans ~= S!bool("--effect"       ,"Perform effect scan", effect);
  opt.booleans ~= S!bool("--qtl"          ,"Perform QTL scan", qtl);
  opt.booleans ~= S!bool("--sp"           ,"Single permutation mode (reuse permutations of T0)", sp);
  opt.strings  ~= S!string("--output"     ,"Path to write output to (DEFAULT: ./)", output);
  opt.booleans ~= S!bool("--redo"         ,"Overwrite previous output files", redo);
  opt.integers ~= S!uint("--nperms"       ,"Number of permutations", nperms);
  opt.doubles  ~= S!double("--minlod"     ,"Minimum LOD score", minlod);
  opt.strings  ~= S!string("--phenotypes" ,"File containing phenotypes", phenotype_filename);
  opt.strings  ~= S!string("--genotypes"  ,"File containing genotypes", genotype_filename);
  opt.strings  ~= S!string("--format"     ,"File format", file_format);
  settings.load(opt);
  return settings;
}
