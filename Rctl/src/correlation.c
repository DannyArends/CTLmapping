/******************************************************************//**
 * \file Rctl/src/correlation.c
 * \brief Implementation of functions related to correlation
 *
 * <i>Copyright (c) 2010-2013</i> GBIC - Danny Arends<br>
 * Last modified Feb, 2013<br>
 * First written 2011<br>
 **********************************************************************/
#include "correlation.h"

void R_correlation(double* x, double* y, double* res, int* dim, int* verb){
  size_t dimension  = (size_t)(*dim);
  bool verbose   = (bool)(*verb);
  res[0] = correlation(x, y, dimension, verbose);
}

void R_correlation1toN(double* x, double* y, double* res, int* dim, int* numy, int* verb){
  size_t i = 0;
  size_t dimension  = (size_t)(*dim);
  size_t ny    = (size_t)(*numy);
  bool verbose = (bool)(*verb);
  double** ynew = asdmatrix(ny, dimension, y);
  double*  cors = cor1toN(x, ynew, dimension, ny, verbose);
  for(i = 0; i < ny; i++){ res[i] = cors[i]; }
  free(cors);
}

void R_chiSQN(int* nr, double* r, double* res, int* phe, int* nsamples, int* nphe){
  size_t p = 0;
  size_t phenotype = (size_t)(*phe);
  size_t nphenotypes = (size_t)(*nphe);
  size_t ncorrelations = (size_t)(*nr);
  double** correlations = asdmatrix(ncorrelations, nphenotypes, r);

  double* chisq = chiSQN(ncorrelations, correlations, phenotype, nsamples, nphenotypes);
  for(p = 0; p < nphenotypes; p++){ if(phenotype != p){ res[p] = chisq[p]; } }
  free(chisq);
}

void R_chiSQtoP(double* Cv, int* Dof, double* res){ res[0] = chiSQtoP((double)(*Cv), (int)(*Dof)); }

double correlation(const double* x, const double* y, size_t dim, bool verbose){
  size_t i;
  double XiYi = 0.0, Xi = 0.0, Yi = 0.0, XiP2 = 0.0, YiP2 = 0.0;
  double onedivn = 1.0 / dim;

  for(i = 0; i < dim; i++){
    if(x[i] != MISSING && y[i] != MISSING){
      XiYi += x[i] * y[i];
      Xi   += x[i];
      Yi   += y[i];
      XiP2 += x[i] * x[i];
      YiP2 += y[i] * y[i];
    }
  }
  double nom = (XiYi - (onedivn *Xi * Yi));
  double denom = sqrt(XiP2 - onedivn * pow(Xi, 2.0)) * sqrt(YiP2 - onedivn * pow(Yi, 2.0));
  double cor = nom / denom;
  if(isNaN(cor) || isinf(cor) || cor < (-1.0 - EPSILON) || cor > (1.0 + EPSILON)){ 
    err("Correlation '%.4f' not in range [-1, 1]\n", cor); 
  }
  return(cor);
}

double* cor1toN(double* x, double** y, size_t dim, size_t ny, bool verbose){
  size_t i, j;
  double nom, denom;
  double onedivn = (1.0 / dim), Xi = 0.0, XiP2 = 0.0;

  double* cors   = newdvector(ny);
  double* Yi     = newdvector(ny);
  double* YiP2   = newdvector(ny);
  double* XiYi   = newdvector(ny);

  // Unrolled 1:N correlation loop
  for(j = 0; j < ny; j++){   // Loop over all traits
    for(i = 0; i < dim; i++){ if(y[j][i] != MISSING && x[i] != MISSING){ // If both are available
      if(j==0){
        Xi   += x[i];
        XiP2 += x[i] * x[i];
      }
      XiYi[j] += x[i] * y[j][i];
      Yi[j]   += y[j][i];
      YiP2[j] += y[j][i] * y[j][i];
    }}
  }

  for(j = 0; j < ny; j++){
    nom   = (XiYi[j] - (onedivn*Xi*Yi[j]));
    denom = sqrt(XiP2 - (onedivn * Xi * Xi)) * sqrt(YiP2[j] - (onedivn * Yi[j] * Yi[j]));
    if(denom == 0) err("Denominator = 0 in correlation (Too few samples in a genotype)\n");
    cors[j] = nom / denom;
    if(isNaN(cors[j]) || isinf(cors[j]) || cors[j] < (-1.0 + EPSILON) || cors[j] > (1.0 + EPSILON)){ 
      err("Correlation '%.8f' not in range [-1, 1]\n", cors[j]);
    }
  }
  free(Yi); free(YiP2); free(XiYi);
  return(cors);
}

double* getCorrelations(const Phenotypes phenotypes, const Genotypes genotypes, size_t phe1, 
                        clvector genoenc, size_t mar, size_t phe2, bool verbose){

  size_t  i;
  double* cors = newdvector(genoenc.nelements);
  if(phe1 != phe2){
    for(i = 0; i < genoenc.nelements; i++){
      clvector inds = which(genotypes.data[mar], phenotypes.nindividuals, genoenc.data[i]);
      double* P1  = get(phenotypes.data[phe1], inds);
      double* P2  = get(phenotypes.data[phe2], inds);
      cors[i]    = correlation(P1, P2, inds.nelements, false);
      if(verbose){
        info("CTL: %d %d %d | %d [%d] -> %f\n", phe1, mar, phe2, genoenc.data[i], inds.nelements, cors[i]);
      }
      free(P1), free(P2); // Clear phenotypes
      free(inds.data);    // Clear index data
      updateR(0);
    }
  }
  return cors;
}

double* chiSQN(size_t nr, double** r, size_t phe, int* nsamples, size_t nphe){
  size_t p, i, denom;
  double sumOfSquares = 0.0, squaresOfSum = 0.0, df;
  double* ret = newdvector(nphe);  /*!< Returned Chi^2 values for phenotype phe against the other phenotypes */
  for(p = 0; p < nphe; p++){
    if(phe != p){
      denom = 0; sumOfSquares = 0.0; squaresOfSum = 0.0; // Reset for next calculation
      for(i = 0; i < nr; i++){
        df = nsamples[i]-3;
        if(df > 0){
          sumOfSquares += df * pow(zscore(r[i][p]), 2.0);
          squaresOfSum += df * zscore(r[i][p]);
          denom += df;
        }
      }
      if(denom == 0) err("Divide by 0 groups too small");
      ret[p] = sumOfSquares - (pow(squaresOfSum, 2.0) / denom);
    }
  }
  return ret;
}

double chiSQtoP(double Cv, int Dof){
  if(Cv <= 0 || Dof < 1) return 1.0;
  return dchisq(Cv,(double)Dof, 0);
}

