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
  void load(){
    settings.booleans ~= S!bool("--help","Shows the help file", false);
    settings.booleans ~= S!bool("--verbose","Verbose mode", true);
    settings.integers ~= S!uint("--nperms","Number of permutations", 100);
    settings.strings  ~= S!string("--pheno","Input filename", "phenotypes.csv");
    settings.strings  ~= S!string("--geno","Input filename", "genotypes.csv");  
  }
  Scontainer settings;
}

CTLsettings parseCmd(string[] args){
  CTLsettings settings;
  settings.load();
  if(args.length > 1){
    foreach(string arg; args[1..$]){
      bool known = false;
      if(arg.indexOf("=") > 1){
        foreach(int cnt, S!string s; settings.settings.strings){
          if(arg[1]==s.name[2]){
            known = true;
            writeln(s.description);
          }
        }
        foreach(int cnt, S!uint s; settings.settings.integers){
          if(arg[1]==s.name[2]){
            known = true;
            writeln(s.description);
          }
        }
      }else{ // Boolean perhaps ?
        foreach(int cnt, S!bool s; settings.settings.booleans){
          if(arg[1]==s.name[2]){
            known = true;
            writeln(s.description);
            settings.settings.booleans[cnt].value = !settings.settings.booleans[cnt].value;
          }
        }
      }
      if(!known) writeln("Unknown argument: ",arg);
    }
  }
  return settings;  
}
