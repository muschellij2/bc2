#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
void getXY( double *x,double *pxx,double *ret ,int *pN,int *pM){
	int NM;
	int N = *pN;
	int M = *pM;
	NM = M*N;
	for(int i=0;i< M;i++){
		ret[i] =0;
		for(int k=0;k<N;k++){
			ret[i] += x[k]*pxx[N*i+k];
		}
	}

return;
}
