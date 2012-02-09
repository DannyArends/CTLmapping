/**********************************************************************
 * src/ctl/core/qtl/nrc.d
 *
 * copyright (c) 2012 Danny Arends
 * last modified Jan, 2012
 * first written Jan, 2012
 **********************************************************************/
module ctl.core.qtl.nrc;

import std.stdio;
import std.math;

/* Standalone functions gammln, betacf, betai necessary to calculate 
   F(P,df1,df2) from Numerical Recipes in C 
   In theory we can get rid of all these function by binding versus R
 */
double gammln(double xx){
  double x,tmp,ser;
  static double cof[6]=[76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5];
  x=xx-1.0;
  tmp=x+5.5;
  tmp-= (x+0.5)*log(tmp);
  ser=1.0;
  for(int j=0; j<=5; j++){ x+=1.0; ser += cof[j]/x; }
  return -tmp+log(2.50662827465*ser);
}

double betacf(double a, double b, double x){
  double qap,qam,qab,em,tem,d,bz,bm=1.0,bp,bpp,az=1.0,am=1.0,ap,app,aold;
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  bz=1.0-qab*x/qap;
  for(int m=1; m<=100; m++){
    em=cast(double)m;
    tem=em+em;
    d=em*(b-em)*x/((qam+tem)*(a+tem));
    ap=az+d*am;
    bp=bz+d*bm;
    d= -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
    app=ap+d*az;
    bpp=bp+d*bz;
    aold=az;
    am=ap/bpp;
    bm=bp/bpp;
    az=app/bpp;
    bz=1.0;
    if(abs((az-aold)/az)  < 3.0e-7) return az;
  }
  writeln("ERROR: a or b too big or max number of iterations too small");
  assert(0);
}

double betai(double a, double b, double x){
  double bt;
  if(x<0.0 || x>1.0){ 
    writeln("ERROR: x not between 0 and 1: ",x);
    assert(0);
  }else{
    if (x==0.0 || x==1.0) bt=0.0;
    else bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
    if (x<(a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
    else return 1.0-bt*betacf(b,a,1.0-x)/b;
  }
}

double inverseF(int df1, int df2, double alfa, bool verbose = false){
  double prob=0.0, minF=0.0, maxF=100.0, halfway=50.0, absdiff=1.0;
 int count=0;
  while ((absdiff>0.001)&&(count<100)){
    count++;
    halfway= (maxF+minF)/2.0;
    prob= betai(df2/2.0,df1/2.0,df2/(df2+df1*halfway));
    if (prob<alfa) maxF= halfway;
    else minF= halfway;
    absdiff= fabs(prob-alfa);
  }
  if(verbose) writeln("prob=" , prob , "; alfa=" , alfa);
  return halfway;
}
