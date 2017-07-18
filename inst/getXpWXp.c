#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void getXpWXp(double *W,double *xp,int *pN,int *pM,double *ret,int *pp){

	int N = *pN;
	int M = *pM;
	int p = *pp;
	for(int l1=0;l1<p;l1++){
		for(int l2=0;l2<p;l2++){
			for(int i=0;i<M;i++){
				for(int j=0;j<M;j++){
				ret[l1*p*M+l2+i*p*p*M+j*p] = 0;
				for(int k=0;k<N;k++){
					ret[l1*p*M+l2+i*p*p*M+j*p] += xp[l1*N+k]*W[(i*M+j)*N+k]*xp[l2*N+k];
				}
			}
		}
		}
	}


return;
}
