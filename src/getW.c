#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
void Weighted_W(double *p, double *W,int *pN,int *pM){
	int NM;
	int N = *pN;
	int M = *pM;
	NM = M*N;
	for(int i=0;i<M;i++){
		for(int j=0;j<M;j++){
		if(i==j){
			for(int k=0;k<N;k++){
				W[NM*i+N*j+k] = p[N*i+k]-p[N*i+k]*p[N*i+k];
			}
		}else{
			for(int k=0;k<N;k++){
				W[NM*i+N*j+k] = -p[N*i+k]*p[N*j+k];
			}
			}

	}
}

return;
}
