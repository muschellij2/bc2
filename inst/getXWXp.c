#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void getXWXp(double *WXp,double *X,int *pN,int *pM,double *ret,int *pp){

	int N = *pN;
	int M = *pM;
	int NM2 = N*M*M;
	int p = *pp;

	for(int l=0;l<p;l++){
		for(int i=0;i<(M*M);i++){
			ret[l*M*M+i] = 0;
			for(int k=0;k<N;k++){
				ret[l*M*M+i] += X[k]*WXp[l*NM2+i*N+k];
			}
		}
		}

return;
}
