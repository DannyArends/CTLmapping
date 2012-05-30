/******************************************************************//**
 * \file ctl/io/qtab/wrapper.d
 * \brief Input reader wrapper for the QTAB format
 *
 * <i>Copyright (c) 2012</i>GBIC - Danny Arends<br>
 * Last modified May, 2012<br>
 * First written Jan, 2012<br>
 * Written in the D Programming Language (http://www.digitalmars.com/d)
 **********************************************************************/
module ctl.io.qtab.wrapper;

version(QTAB){

import std.stdio, std.conv;
import ctl.core.array.matrix, ctl.io.reader, ctl.io.terminal;
import qtl.core.phenotype, qtl.core.genotype, qtl.core.primitives;
import qtl.plugins.qtab.read_qtab;

class QTABreader : Reader{
  double[][] loadphenotypes(string filename = "phenotype.qtab"){
    MSG("Starting with qtab phenotypes");
    auto res = read_phenotype_qtab!(Phenotype!double)(filename);
    MSG("Done with qtab phenotypes");
    Phenotype!double[][] pheno = res[0];
    double[][] phenotypes = newmatrix!double(pheno.length,pheno[0].length, 0.0);
    for(size_t x=0;x< pheno.length;x++){
      for(size_t y=0;y< pheno[0].length;y++){
        if(pheno[x][y].value != double.max) phenotypes[x][y] = pheno[x][y].value;
      }
    }
    return translate(phenotypes);
  }
    
  int[][] loadgenotypes(string filename = "genotype"){
    string symbol_fn   = filename ~ "_symbols.qtab";  // FIXME: file name conventions are funny
    string genotype_fn = filename ~ "_genotypes.qtab";
    MSG("Starting with qtab genotypes");
    MSG("Reading " ~ symbol_fn);
    auto symbols = read_genotype_symbol_qtab(File(symbol_fn,"r"));
    // assert(to!string(symbols.decode("A")) == "[(0,0)]");
    MSG("Done with qtab genotypes symbols");
    MSG("Reading " ~ genotype_fn);
    auto ret = read_genotype_qtab(File(genotype_fn,"r"), symbols);
    MSG("Done with qtab genotypes");
    auto individuals = ret[0];
    auto genotype_matrix = ret[1];
    
    int[][] genotypes = newmatrix!int(genotype_matrix.length,genotype_matrix[0].length, 0);
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
    return translate(genotypes);
  }
}

} // version(QTAB)
