#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void getXWX(double *W,double *X,int *pN,int *pM,double *ret){

	int N = *pN;
	int M = *pM;
	//scint NM = N*M;
	double * Xtest2;
	Xtest2 = (double*) malloc(N*sizeof(double));


	for(int i=0;i<N;i++){
		Xtest2[i] = X[i]*X[i];
	}

	for(int i=0;i<M*M;i++){
		ret[i] =0;
		for(int k=0;k<N;k++){
			ret[i] += Xtest2[k]*W[i*N+k];
		}
	}
	free(Xtest2);
return;
}
