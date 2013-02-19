/******************************************************************//**
 * \file Rctl/src/correlation.c
 * \brief Implementation of functions related to correlation
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "correlation.h"

double correlation(const double* x, const double* y, size_t dim, bool verbose){
  size_t i;
  double onedivn = 1.0 / dim;
  KahanAccumulator XiYi = createAccumulator();
  KahanAccumulator Xi   = createAccumulator();
  KahanAccumulator Yi   = createAccumulator();
  KahanAccumulator XiP2 = createAccumulator();
  KahanAccumulator YiP2 = createAccumulator();

  for(i = 0; i < dim; i++){
    if(x[i] != MISSING && y[i] != MISSING){
      XiYi = KahanSum(XiYi, x[i] * y[i]);
      Xi   = KahanSum(Xi, x[i]);
      Yi   = KahanSum(Yi, y[i]);
      XiP2 = KahanSum(XiP2, x[i] * x[i]);
      YiP2 = KahanSum(YiP2, y[i] * y[i]);
    }
  }
  double nom = (XiYi.sum - (onedivn*Xi.sum*Yi.sum));
  double denom = sqrt(XiP2.sum - onedivn * pow(Xi.sum, 2.0)) * sqrt(YiP2.sum - onedivn * pow(Yi.sum, 2.0));
  double cor = nom/denom;
  if(isNaN(cor) || isinf(cor) || cor < -1 || cor > 1){ 
    err("Correlation '%.4f' not in range [-1, 1]\n", cor); 
  }
  return(cor);
}

double* cor1toN(const double* x, double** y, size_t dim, size_t ny, bool verbose){
  size_t i, j;
  double  onedivn = 1.0 / dim;

  double  Xi   = 0.0;
  double  XiP2 = 0.0;

  double* Yi   = newdvector(ny);
  double* YiP2 = newdvector(ny);
  double* XiYi = newdvector(ny);

  for(i = 0; i < dim; i++){ if(x[i] != MISSING){ // Loop over the non-missing in the common dimension
    Xi   += x[i];
    XiP2 += x[i] * x[i];
    for(j = 0; j < ny; j++){ if(y[j][i] != MISSING){ // Same for Y
      XiYi[j] += x[i] * y[j][i];
      Yi[j]   += y[j][i];
      YiP2[j] += y[j][i] * y[j][i];
    }}
  }}
  double* cors = newdvector(ny);
  for(j = 0; j < ny; j++){
    double nom = (XiYi[j] - (onedivn*Xi*Yi[j]));
    double denom = sqrt(XiP2 - onedivn * pow(Xi, 2.0)) * sqrt(YiP2[j] - onedivn * pow(Yi[j], 2.0));
    cors[j] = 1e-6 * (int)((nom / denom) * 1e6);  // TODO: Remove this easy fix for rounding errors
    if(isNaN(cors[j]) || isinf(cors[j]) || cors[j] < -1.0 || cors[j] > 1.0){ 
      err("Correlation '%.8f' not in range [-1, 1]\n", cors[j]);
    }
  }
  free(Yi); free(YiP2); free(XiYi);
  return(cors);
}

double* getCorrelations(const Phenotypes phenotypes, const Genotypes genotypes, size_t phe1, 
                        clvector genoenc, size_t mar, size_t phe2, bool verbose){

  size_t  i;
  double* cors  = newdvector(genoenc.nelements);
  if(phe1 != phe2){
    for(i = 0; i < genoenc.nelements; i++){
      clvector inds = which(genotypes.data[mar], phenotypes.nindividuals, genoenc.data[i]);
      double* P1  = get(phenotypes.data[phe1], inds);
      double* P2  = get(phenotypes.data[phe2], inds);
      cors[i]    = correlation(P1, P2, inds.nelements, false);
      if(verbose){
        info("Significant: %d %d %d | %d [%d] -> %f\n", phe1, mar, phe2, genoenc.data[i], inds.nelements, cors[i]);
      }
      free(P1), free(P2); // Clear phenotypes
      free(inds.data);    // Clear index data
      updateR(0);
    }
  }
  return cors;
}

double* chiSQN(size_t nr, double** r, size_t phe, int* nsamples, size_t nphe){
  size_t p, i;
  double* ret = newdvector(nphe);
  for(p = 0; p < nphe; p++){
    if(phe != p){
      double* ts = newdvector(nr);
      for(i = 0; i < nr; i++){
        ts[i] = r[i][p];
      }
      ret[p] = chiSQ(nr, ts, nsamples);
      free(ts);
    }
  }
  return ret;
}

double chiSQ(size_t nr, double* r, int* nsamples){
  size_t i;
  size_t denom = 0;
  KahanAccumulator sumOfSquares = createAccumulator();
  KahanAccumulator squaresOfSum = createAccumulator();

  for(i = 0; i < nr; i++){
    sumOfSquares = KahanSum(sumOfSquares, (nsamples[i]-3) * pow(zscore(r[i]), 2.0));
    squaresOfSum = KahanSum(squaresOfSum, (nsamples[i]-3) * zscore(r[i]));
    denom  += (nsamples[i]-3);
  }
  if(denom == 0) err("Divide by 0 groups too small");
  return(sumOfSquares.sum - (pow(squaresOfSum.sum, 2.0) / denom));
}

double chiSQtoP(int Dof, double Cv){
  if(Cv <= 0 || Dof < 1) return 1.0;
  return dchisq(Cv,(double)Dof, 0);
}

