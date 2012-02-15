/**********************************************************************
 * src/ctl/io/qtab/wrapper.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Feb, 2012
 * first written May, 2011
 **********************************************************************/
module ctl.io.qtab.wrapper;

version(QTAB) {

import std.stdio;
import std.conv;

import ctl.core.array.matrix;
import ctl.io.reader;

import qtl.plugins.qtab.read_qtab;
import qtl.core.phenotype;
import qtl.core.genotype;
import qtl.core.primitives;

class QTABreader : Reader{
  double[][] loadphenotypes(string filename = "phenotype.qtab"){
    writeln("Starting with qtab phenotypes");
    auto res = read_phenotype_qtab!(Phenotype!double)(filename);
    writeln("Done with qtab phenotypes");
    Phenotype!double[][] pheno = res[0];
    double[][] phenotypes = newmatrix!double(pheno.length,pheno[0].length);
    for(size_t x=0;x< pheno.length;x++){
      for(size_t y=0;y< pheno[0].length;y++){
        phenotypes[x][y] = pheno[x][y].value;
      }
    }
    return phenotypes;
  }
    
  int[][] loadgenotypes(string filename = "genotype"){
    string symbol_fn   = filename ~ "_symbols.qtab";  // FIXME: file name conventions are funny
    string genotype_fn = filename ~ "_genotypes.qtab";
    writeln("Starting with qtab genotypes");
    writeln("Reading " ~ symbol_fn);
    auto symbols = read_genotype_symbol_qtab(File(symbol_fn,"r"));
    // assert(to!string(symbols.decode("A")) == "[(0,0)]");
    writeln("Done with qtab genotypes symbols");
    writeln("Reading " ~ genotype_fn);
    auto ret = read_genotype_qtab(File(genotype_fn,"r"), symbols);
    writeln("Done with qtab genotypes");
    auto individuals = ret[0];
    auto genotype_matrix = ret[1];
    
    int[][] genotypes = newmatrix!int(genotype_matrix.length,genotype_matrix[0].length);
    
    for(size_t x=0;x< genotype_matrix.length;x++){
      for(size_t y=0;y< genotype_matrix[0].length;y++){
      if(genotype_matrix[x][y] == symbols.decode("A")){
        genotypes[x][y] = 0;
      }
      if(genotype_matrix[x][y] == symbols.decode("B")){
        genotypes[x][y] = 1;
      }
      }
    }
    return genotypes;
  }
}

} // version(QTAB)

