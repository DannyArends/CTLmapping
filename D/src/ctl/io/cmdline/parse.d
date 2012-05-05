/**********************************************************************
 * src/ctl/io/cmdline/parse.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.io.cmdline.parse;

import std.stdio;
import std.string;
import std.conv;

struct S(T){
  string name;
  string description;
  T      value;
}

struct Scontainer{
  S!bool[]   booleans;
  S!uint[]   integers;
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
      foreach(S!string s; opt.strings){ writeln("  " ~ s.name[0] ~ "[" ~ s.name[1..3] ~"]" ~ s.name[3..$] ~ " - " ~ s.description); }
      return true;
    }
    return false;
  }
  
  uint getInt(string name){
    foreach(int cnt, S!uint s; opt.integers){
      if(name==s.name){
        return s.value;
      }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return -1;
  }
  
  bool getBool(string name){
    foreach(int cnt, S!bool s; opt.booleans){
      if(name==s.name){
        return s.value;
      }
    }
    writeln(" - ERROR: Unknown parameter requested: ",name);
    return false;
  }
  
  string getString(string name){
    foreach(int cnt, S!string s; opt.strings){
      if(name==s.name){
        return s.value;
      }
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
  bool verbose = false;
  bool overwrite = false;
  uint nperms = 100;
  string phenotype_filename = "test/data/phenotypes.csv";
  string genotype_filename = "test/data/genotypes.csv";
  string file_format = "csv";

  getopt(args, "help|h", &help
             , "verbose|v", &verbose
             , "overwrite|o", &overwrite
             , "nperms|n", &nperms
             , "phenotypes|p", &phenotype_filename
             , "genotypes|g", &genotype_filename
             , "format|f", &file_format);

  opt.booleans ~= S!bool("--help"         ,"Show the help file", help);
  opt.booleans ~= S!bool("--verbose"      ,"Verbose mode", verbose);
  opt.booleans ~= S!bool("--overwrite"    ,"Overwrite previous output files", overwrite);
  opt.integers ~= S!uint("--nperms"       ,"Number of permutations", nperms);
  opt.strings  ~= S!string("--phenotypes" ,"File containing phenotypes", phenotype_filename);
  opt.strings  ~= S!string("--genotypes"  ,"File containing genotypes", genotype_filename);
  opt.strings  ~= S!string("--format"     ,"File format", file_format);
  settings.load(opt);
  return settings;
}

