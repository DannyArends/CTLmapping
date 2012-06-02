/******************************************************************//**
 * \file ctl/mapctl.d
 * \brief Function holding the main mapctl workflow
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written May, 2011<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
import std.stdio, std.math, std.conv, std.file, std.datetime;
import ctl.core.array.matrix, core.memory;
import ctl.core.stats.basic, ctl.core.stats.tolod;
import ctl.core.analysis, ctl.core.ctl.utils;
import ctl.core.ctl.mapping, ctl.core.ctl.permutation;
import ctl.io.reader, ctl.io.terminal;
import ctl.io.cmdline.parse;
import ctl.io.csv.write, ctl.io.csv.parse;

void main(string[] args){
  SysTime stime = Clock.currTime();
  MSG("Correlated Trait Locus (CTL) mapping in D");
  MSG("(c) 2012 written by Danny Arends in the D programming language");
  CTLsettings settings   = parseCmd(args);
  MSG("ARGS parsed");
  Reader ireader  = initialize(settings);
  bool verbose    = settings.getBool("--verbose");
  bool overwrite  = settings.getBool("--redo");
  string output   = settings.getString("--output");
  if(output[($-1)]=='/') output = output[0..($-1)];
  string input_p  = settings.getString("--phenotypes");
  string input_g  = settings.getString("--genotypes");

  if(overwrite) WARN("Overwriting files in on");
  if(verbose){ MSG("Verbose mode on"); MSG("Output saved to: " ~ output); }
  if(!settings.displayHelp()){
    MSG("Start loading input files (%s, %s)",  input_p, input_g);
    double[][]  phenotypes = ireader.loadphenotypes(input_p);
    string[]    phenonames = ireader.loadphenonames(input_p);
    int[][]     genotypes  = ireader.loadgenotypes(input_g);
    if(verbose){
      MSG("Dataset: %s geno- and %s phenotypes", genotypes.length, phenotypes.length);
      if(genotypes.length == 0){ ERR("No genotypes loaded, analysis aborted"); return; }
      if(phenotypes.length == 0){ ERR("No phenotypes loaded, analysis aborted"); return; }
      MSG("Dataset: Measurements on %s and %s individuals\n", genotypes[0].length, phenotypes[0].length);
    }
    if(genotypes[0].length != phenotypes[0].length){ ERR("Mismatch between individuals %s != %s", genotypes[0].length, phenotypes[0].length); return; }
    //We need an output path default to ./ ??
    if(!exists(output)) mkdirRecurse(output);

    double[][] effects, qtllod, ctllod, score, perms;

    Analysis analysis;
    //Start by mapping all QTL
    if(needanalysis(output ~ "/effects.txt",overwrite)){
      analysis = getanalysis("effect",settings);
      effects = analysis.analyse(genotypes, phenotypes, [], verbose);
      if(effects != null) writeFile(effects, output ~ "/effects.txt", null, overwrite, verbose);
    }else{ WARN("Skipped effect scan"); }

    //Start by mapping all QTL
    if(needanalysis(output ~ "/qtls.txt",overwrite)){
      analysis = getanalysis("qtl",settings);
      qtllod = analysis.analyse(genotypes, phenotypes, [], verbose);
      if(qtllod != null) writeFile(qtllod, output ~ "/qtls.txt", null, overwrite, verbose);
    }else{ WARN("Skipped QTL mapping"); }

    for(size_t p=0; p < phenotypes.length; p++){ //Main CTL mapping loop
      if(verbose) MSG("- Phenotype %s -",p);
      string fn_ctl  = output ~ "/ctl"~to!string(p)~".txt";
      string fn_perm = output ~ "/perms"~to!string(p)~".txt";
      string fn_lods = output ~ "/lodscores"~to!string(p)~".txt";

      if(needanalysis(fn_ctl,overwrite)){
        score = mapping(phenotypes,  genotypes, p, true);
        writeFile(translate(score),  fn_ctl, null, overwrite, verbose);
      }else{ //Reload the scores
        score = parseFile!double(fn_ctl, false, false);
        MSG("Skipped CTL mapping, file %s exists", fn_ctl); }
     if(settings.getBool("--sp") && p > 0){
        MSG("Reusing permutations");
      }else{
        if(needanalysis(fn_perm,overwrite)){
          perms = permutation(phenotypes, genotypes, p, settings.getInt("--nperms"), verbose);
          writeFile(translate(perms), fn_perm, null, overwrite, verbose);
        }else{ //Reload the permutations
          perms = parseFile!double(fn_ctl, false, false);
          MSG("Skipped permutations, file %s exists", fn_perm); }
      }
      if(needanalysis(fn_lods,overwrite)){
        double[] permlist    = doMatrixMax(absmatrix(perms));
        MSG("Permlist: %s", permlist.length);
        score                = translate(score);
        permlist.sort;
        double   cutoff      = getCutoff(permlist, settings.getDouble("--minlod"));
        size_t[] significant = getSignificant(score, cutoff);
        ctllod = tolod(score, permlist, significant, phenonames, verbose);
        writeFile(ctllod, fn_lods, get(phenonames,significant), overwrite, verbose);
      }else{ MSG("Skipped LOD transformation, file %s exists", fn_lods); }
      GC.collect();
      GC.minimize();
    }
    writeln();
    MSG("mapCTL finished (%s seconds)",(Clock.currTime()-stime).total!"seconds"()," seconds");
    MSG("Continue by starting R and loading the results:\n library(ctl)\n ctls <- load.ctl(\"%s\", \"%s\", \"%s\")\n plot(ctls)",input_g, input_p,output);
  }
}
