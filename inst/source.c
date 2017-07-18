#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define ALMOST_ZERO 1e-16
#define NUM_ZERO 1e-100
#define ERROR_SINGULAR_MATRIX 1
#define CHECK_MEM(obj) if (obj == NULL) {Rprintf("ERROR: allocating memory \n"); error("1");}

static void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
static void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %d ", vec[i]);
  }
  printf("\n \n");
}


static void print_dMat(mat, nr, nc, name)
double **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %g ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}

static void print_iMat(mat, nr, nc, name)
int **mat;
int nr, nc;
char name[10];
{
  int i, j;
  printf("%s \n", name);
  for (i=0; i<nr; i++) {
    for (j=0; j<nc; j++) printf(" %d ", mat[i][j]);
    printf("\n");
  }
  printf("\n \n");
}

/* Function to allocate memory for a double vector */
static double * dVec_alloc(n, initFlag, initVal)
int n, initFlag;
double initVal;
{
  int i;
  double *ret, *p;

  ret = (double *) malloc(n*sizeof(double));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }
  
  return(ret);

} /* END: dVec_alloc */

/* Function to allocate a double matrix */
static double ** dMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
double initVal;
{
  double **mat, **ptr;
  int i;

  mat = (double **) malloc(nrow*sizeof(double *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = dVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: dMat_alloc */

/* Function to allocate memory for a integer vector */
static int * iVec_alloc(n, initFlag, initVal)
int n, initFlag;
int initVal;
{
  int i;
  int *ret, *p;

  ret = (int *) malloc(n*sizeof(int));
  CHECK_MEM(ret);
  if (initFlag) {
    for (i=0, p=ret; i<n; i++, p++) *p = initVal;
  }
  
  return(ret);

} /* END: iVec_alloc */

/* Function to allocate a int matrix */
static int ** iMat_alloc(nrow, ncol, initFlag, initVal)
int nrow, ncol, initFlag;
int initVal;
{
  int **mat, **ptr;
  int i;

  mat = (int **) malloc(nrow*sizeof(int *));
  CHECK_MEM(mat);
  for (i=0, ptr=mat; i<nrow; i++, ptr++) *ptr = iVec_alloc(ncol, initFlag, initVal);

  return(mat);

} /* END: iMat_alloc */

/* Function to free a matrix */
static void matrix_free(x, n)
void **x;
int n;
{
  int i;
  for (i=0; i<n; i++) {
    if (x[i]) free(x[i]);
  }
  free(x);

} /* END: matrix_free */

/* Function to check for non-finite values */
static int all_finite(vec, n)
double *vec;
int n;
{
  int i;

  for (i=0; i<n; i++) {
    if (!R_FINITE(vec[i])) return(0);
  } 

  return(1);

} /* END: all_finite */

/* Multiply to matrices (matrix mult)  */
static void matrixMult(m1, m1_nr, m1_nc, m2, m2_nc, ret)
double **m1, **m2, **ret;
int m1_nr, m1_nc, m2_nc;
{
  int i, j, k;
  double sum;

  for (i=0; i<m1_nr; i++) {
    for (j=0; j<m2_nc; j++) {
      sum = 0.;
      for (k=0; k<m1_nc; k++) sum += m1[i][k]*m2[k][j]; 
      ret[i][j] = sum;
    }
  }

} /* END: matrixMult */

/* Function to compute Xy for matrix X and vector */
static void X_y(X, nr, nc, y, ret)
double **X, *y, *ret;
int nr, nc;
{
  int i, j;
  double sum, *p, *pret, *px;

  for (i=0, pret=ret; i<nr; i++, pret++) {
    sum  = 0.0;
    for (j=0, p=y, px=X[i]; j<nc; j++, p++, px++) {
      sum += *px * *p;
    }
    *pret = sum; 
  }

} /* END: X_y */

/* Function for dot product of two vectors */
static double dotProd(v1, v2, n)
double *v1, *v2;
int n;
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += v1[i]*v2[i];

  return(sum);

} /* END: dotProd */

/* Function to compute Z */

static void get_Z(Z_design, M, Ncatp1, Ncov, out)
double **Z_design, **out;
int M, Ncatp1, Ncov;
{
  int i, j, k, row=0, col;
  double *vec;


  /* Get the column for diag(P) (X) z_design, loop over rows of z_design */
  for(i = 0;i<M;i++){
      row = 0+ i*Ncov;
      col = i;
      out[row][col] = 1.0;
  }
  for(i = 0;i<M;i++){
    vec = Z_design[i];
    for(j =0; j<(Ncov-1);j++){
       row = j+1+i*(Ncov);
       col = M+j*Ncatp1;
      for(k=0;k<Ncatp1;k++){
      out[row][col] = vec[k];
      col++;
    }
    }

  }

} /* END: get_Z */

/* Function to get a column of XXZ matrix */
static void get_XXZ_col(col, X, Xnr, Xnc, M, Z, dtempvec, out)
double **X, **Z, *out, *dtempvec; /* dtempvec length M*Xnc */
int col, Xnr, Xnc, M;
{
  int  i, j, k, zstart;
  double sum, *pt, *prow, *pout;

  /* Copy column of Z into dtempvec */
  for (i=0, pt=dtempvec; i<Xnc*M; i++, pt++) *pt = Z[i][col];

  pout = out;
  for (i=0; i<M; i++) {
    zstart = i*Xnc;
    for (j=0; j<Xnr; j++) {
      sum  = 0.0;
      for (k=0, prow=X[j], pt=&dtempvec[zstart]; k<Xnc; k++, prow++, pt++) {
        sum += *prow * *pt;
      }
      *pout++ = sum;
    }
  }

} /* END: get_XXZ_col */

/* Function to compute t(XXZ) */
static void get_tXXZ(X, Xnr, Xnc, M, Ncatp1, Z, out)
double **X, **Z, **out;  /* out must have dim Xnc*Ncatp1 x Xnr*M */
int Xnr, Xnc, M, Ncatp1;
{
  int i;
  double *temp;

  temp = dVec_alloc(Xnc*M, 0, 0.0);

  for (i=0; i<((Xnc-1)*Ncatp1+M); i++) {
    get_XXZ_col(i, X, Xnr, Xnc, M, Z, temp, out[i]);
  }

  free(temp);

} /* END: get_tXXZ */

/* Function to compute lxx = xx * beta = xx * z * delta */
static void get_lxx(Z, Znr, Znc, delta, X, Xnr, Xnc, M, out)
double **Z, *delta, **X, *out;
int Xnr, Xnc, Znr, Znc, M;
{
  double *beta, sum;
  int i, j, k, row, brow, bstart;

  beta  = dVec_alloc(Znr, 0, 0.0);

  /* Compute beta = Z*delta */ 
  X_y(Z, Znr, Znc, delta, beta);  

  /* Compute out = XX*beta */
  row = 0;
  for (i=0; i<M; i++) {
    bstart = i*Xnc;
    for (j=0; j<Xnr; j++) {
      sum  = 0.0;
      brow = bstart;
      for (k=0; k<Xnc; k++) {
        sum += X[j][k]*beta[brow];
        brow++;
      }
      out[row] = sum;
      row++;
    }
  }

  free(beta);

} /* END: get_lxx */

/* Function to compute pxx */
static void get_pxx(lxx, N, M, out)
double *lxx, *out;
int N, M;
{
  int i, j, col, NM;
  double sum;

  NM = N*M;
  for (i=0; i<NM; i++) out[i] = exp(lxx[i]);
   
  /* Scale */
  for (i=0; i<N; i++) {
    sum = 0.0;
    col = i;
    for (j=0; j<M; j++) {
      sum += out[col];
      col += N; 
    }
    col = i;
    sum = sum + 1.0;
    for (j=0; j<M; j++) {
      out[col] = out[col]/sum;
      col += N; 
    }
  }

} /* END: get_pxx */

/* Function to compute vec1*W*vec2 */
static double v1Wv2(p, N, M, vec1, vec2)
double *p, *vec1, *vec2;  /* p is stored as a vector, out must be of length NM */
int N, M;
{
  int i, ii, jj, NM, MP1, row, NMP1;
  double sum, prow, *p1, *pv2, *pv1, ret;

  NM   = N*M;
  MP1  = M + 1;
  NMP1 = NM + 1;
  
  ret = 0.0;
  for (row=0, p1=p, pv2=vec2, pv1=vec1; row<NM; row++, p1++, pv2++, pv1++) {
    prow = *p1;
    sum  = (prow-prow*prow)* *pv2;
    ii   = row + N;
    jj   = row - N;
    for (i=2; i<MP1; i++) {
      if (ii < NMP1) {
        sum += -prow*p[ii]*vec2[ii];
        ii   = ii + N;
      }
      if (jj > -1) {
        sum += -prow*p[jj]*vec2[jj];
        jj   = jj - N;
      }
    }
    ret += *pv1 * sum;
  }

  return(ret);

} /* END: v1Wv2 */

/* Function to compute W_y = yy-pxx+W%*%lxx */
static void get_Wy2(Y, lxx, Pxx, W, W_index, N, M, out)
double *Y;  /* Nsub*M vector of outcomes */
int N, M, **W_index;
double *lxx, *out, *Pxx, **W;
{
  int i, NM, row, **prWi, *pcWi, cnt;
  double sum, *p2, **pr, *pc, *p1, *py;

  NM   = N*M;
  
  prWi = W_index;
  cnt  = 1;
  for (row=0, pr=W, p2=out, py=Y, p1=Pxx; row<NM; row++, pr++, p2++, py++, p1++) {
    sum = 0.0;
    for (i=0, pc=*pr, pcWi=*prWi; i<M; i++, pc++, pcWi++) {    
      sum += *pc * lxx[*pcWi];
    }
    *p2 = *py - *p1 + sum;

    cnt++;
    if (cnt > N) {
      cnt  = 1;
      prWi = W_index;
    } else {
      prWi++;
    }
  }

} /* END: get_Wy2 */

/* Function to compute information matrix t(xxz)%*%W%*%xxz */
static void get_Info2(tXXZ, N, M, nparm, W, W_index, out)
double **tXXZ, **out, **W;
int N, M, nparm, **W_index;
{
  int i, j, row, col, NM, **pr, *pc, cnt;
  double sum, **prW, *pcW, sum2, **pro, *pco, **prx, *pcx, *pcx2;

  NM = N*M;
  pr = W_index;

  for (row=0, pro=out, prx=tXXZ; row<nparm; row++, pro++, prx++) {
    pco = *pro;
    for (col=row; col<nparm; col++) {
      /* Compute XXZ[, row]*W*XXZ[, col] */
      sum2 = 0.0;
      cnt  = 1;
      pcx2 = tXXZ[col];
      for (i=0, prW=W, pcx=*prx; i<NM; i++, prW++, pcx++) {
        sum = 0.0;
        for (j=0, pc=*pr, pcW=*prW; j<M; j++, pc++, pcW++) {
          sum += *pcW * pcx2[*pc];
        } 
        sum2 += *pcx * sum;
        cnt++;
        if (cnt > N) {
          cnt = 1;
          pr  = W_index;
        } else {
          pr++;
        }
      }
      pco[col] = sum2;
      if (row != col) out[col][row] = sum2; 
    }
  }

} /* END: get_Info2 */


/* Function to compute delta */
static void get_delta(INV, nc, tXXZ, N, M, W_y, out)
double **INV, **tXXZ, *W_y, *out; 
int nc, N, M;
{
  /* delta <- solve(Infor_M,t(xxz)%*%W_y) */
  int i, NM;
  double *p1, **p2, *tXXZWy;

  NM     = N*M;
  tXXZWy = dVec_alloc(nc, 0, 0.0);

  /* Compute t(xxz)%*%w_y */
  for (i=0, p1=tXXZWy; i<nc; i++, p1++) {
    *p1 = dotProd(tXXZ[i], W_y, NM);
  }

  /* Compute delta */
  for (i=0, p1=out, p2=INV; i<nc; i++, p1++, p2++) *p1 = dotProd(*p2, tXXZWy, nc);
   
  free(tXXZWy);

} /* END: get_delta */


/* Function to check the stopping criteria */
static double checkStop(delta, delta0, n)
double *delta, *delta0;
int n;
{
  int i;
  double maxv=-9999.9, temp, rerror;

  /* Get max value */
  for (i=0; i<n; i++) {
    temp = fabs(delta0[i]);
    if (temp > maxv) maxv = temp;
  }

  if (maxv < 0.1) maxv = 0.1;

  rerror = -9999.9;
  for (i=0; i<n; i++) {
    temp = fabs(delta[i] - delta0[i])/maxv;
    if (temp > rerror) rerror = temp;
  }

  return(rerror);

} /* END: checkStop */

/* Function to fill in a matrix from a vector (by column) */
static void fillMat(vec, nr, nc, addInt, out)
double *vec, **out;
int nr, nc, addInt;
{
  int i, j, col=0, ii;

  if (addInt) {
    /* Intercept for first column */
    for (i=0; i<nr; i++) out[i][0] = 1.0;
    col = 1;
  }

  ii = 0;
  for (j=0; j<nc; j++) {
    for (i=0; i<nr; i++) {
      out[i][col] = vec[ii];
      ii++; 
    }
    col++;
  }

} /* END: fillMat */

/***********************************************************************************/
/* For computing inverse */
/***********************************************************************************/

/* Factor a symmetric positive definite matrix */
static int symPosFactor(mat, n, ret, retdiag)
double **mat, **ret, *retdiag;
int n;
{
  int i, j, k;
  double sum, save, *ptr;
  
  /* Copy mat to ret */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) ret[i][j] = mat[i][j];
  }

  for (i=0, ptr=retdiag; i<n; i++, ptr++) {
    for (j=i; j<n; j++) {
      sum = ret[i][j];
      for (k=i-1; k>-1; k--) sum -= ret[i][k]*ret[j][k];
      if (i == j) {
        if (sum < NUM_ZERO) return(ERROR_SINGULAR_MATRIX); 
        save = sqrt(sum);
        *ptr = save;
      } else {
        ret[j][i] = sum/save;
      }
    }
  }

  /* Zero out the diagonal and above */
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) ret[i][j] = 0.;
  }

  return(0);

} /* END: symPosFactor */

/* Invert a factor */
static void symPosFacInv(L, diag, n, ret)
double **L, *diag;
int n;
double **ret;
{
  int i, j, k;
  double sum;
  
  /* Copy L to ret */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) ret[i][j] = L[i][j];
  }

  for (i=0; i<n; i++) {
    ret[i][i] = 1./diag[i];
    for (j=i+1; j<n; j++) {
      sum = 0.;
      for (k=i; k<=j-1; k++) sum -= ret[j][k]*ret[k][i];
   
      ret[j][i] = sum/diag[j];
    }
  }

} /* END: symPosFacInv */

/* Transpose of a square matrix */
static void matTranspose(mat, n, ret)
double **mat;
int n;
double **ret;
{
  int i, j;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) ret[j][i] = mat[i][j];   
  }

} /* END: rmatTranspose */

/* Inverse of symmetric positive definite matrix */
static int symPosMatInv(mat, n, ret)
double **mat; 
int n;
double **ret;
{
  double **L, **Linv, *diag;  
  int i;

  L = dMat_alloc(n, n, 1, 0.0);
  diag = (double *) dVec_alloc(n, 1, 0.0);

  /* Get cholesky lower triangle */
  i = symPosFactor(mat, n, L, diag);
  if (i) {
    free(diag);
    matrix_free((void **) L, n);
    return(i);
  }

  /* Get the lower inverse */
  Linv = dMat_alloc(n, n, 1, 0.0);
  symPosFacInv(L, diag, n, Linv);
  free(diag);

  /* Get the transpose of Linv */
  matTranspose(Linv, n, L);

  /* Inverse is t(L^-1)(L^-1) */
  matrixMult(L, n, n, Linv, n, ret);

  matrix_free((void **) L, n);
  matrix_free((void **) Linv, n);

  return(0);

} /* END: symPosMatInv */

/* Function to compute the inverse of a covariance matrix */
int cov_inv(cov, n, inv)
double **cov;
int n;
double **inv; /* Returned inverse */
{
  double cc, a, b, d;
  int ret;

  switch (n) {
    case 0:
      Rprintf("\nERROR: dimension of covariance matrix is 0\n");
      exit(1);
    case 1:
      a = cov[0][0];
      if (fabs(a) < ALMOST_ZERO) return(ERROR_SINGULAR_MATRIX);
      inv[0][0] = 1.0/a; 
      break;
    case 2:
      a  = cov[0][0];
      b  = cov[0][1];
      d  = cov[1][1];
      cc = a*d - b*b;
      if (fabs(cc) < ALMOST_ZERO) return(ERROR_SINGULAR_MATRIX);
      cc = 1.0/cc;
      inv[0][0] = d*cc;
      inv[0][1] = -b*cc;
      inv[1][0] = -b*cc;
      inv[1][1] = a*cc;
      break;
    default:
      ret = symPosMatInv(cov, n, inv);
      if (ret) return(ret);
      break; 
  } /* END: switch */

  return(0);

} /* END: cov_inv */

/* Function to get indices for W matrix */
static void get_W_index(N, M, out)
int N, M, **out; /* out should be NxM */
{
  int j, k, **prow, *pvec, add;

  for (j=0, prow=out; j<N; j++, prow++) {
    add = 0;
    for (k=0, pvec=*prow; k<M; k++, pvec++) {
      *pvec = j + add;
      add   = add + N;
    }
  }

} /* END: get_W_index */

/* Function to compute W */
static void get_W(Pxx, N, M, out)
double *Pxx, **out;
int N, M;
{
  int i, j, k, Ni, Nj, m, n;
  double p;  

  for (i=0; i<M; i++) {
    Ni = N*i;
    for(j=0; j<M; j++) {
       if (i == j) {
         for(k=0; k<N; k++) {
           m = Ni + k;
           p = Pxx[m];
           out[m][i] = p-p*p;
         }
       } else{
         Nj = N*j;
         for(k=0; k<N; k++){
            m = Ni + k;
            n = Nj + k;
            out[m][j] = -Pxx[m]*Pxx[n];
         }
       }
     }
   }

} /* END: get_W */




/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/

/* Function Mvpoly */
static void mvpoly2(delta0, Nparm, Y, X, Z, N, M, Ncat, Ncov, Niter, tol, 
            ret_delta, W, W_index, Inv, pxx, tXXZ, DEBUG)
double *delta0, *Y, **X, **Z, tol, **W, **Inv, *pxx, *ret_delta, **tXXZ;
int Nparm, N, M, Ncat, Ncov, Niter, DEBUG, **W_index;
{
  int i, Znr, Znc, NM, rc, conv=0, Ncatp1, iter;
  double rerror, *w_y, **Info, *lxx;
  
  Ncatp1  = Ncat + 1;
  Znr     = M*Ncov;
  Znc     = (Ncov-1)*Ncatp1+M;
  NM      = N*M;

  /* Allocate memory for matrix of covariates and Z map matrix */
  if (DEBUG) Rprintf("Allocate memory\n");
  lxx      = dVec_alloc(NM, 0, 0.0);
  Info     = dMat_alloc(Nparm, Nparm, 0, 0.0);
  w_y      = dVec_alloc(NM, 0, 0.0);
 
  if (DEBUG) Rprintf("Begin the main loop\n");
  for (iter=1; iter<=Niter; iter++) {
    if (DEBUG) Rprintf("Iteration: %d\n", iter);

    if (DEBUG) Rprintf("Compute lxx\n");
    get_lxx(Z, Znr, Znc, delta0, X, N, Ncov, M, lxx);

    if (DEBUG) Rprintf("Compute pxx\n");
    get_pxx(lxx, N, M, pxx);

    if (DEBUG) Rprintf("Compute W\n");
    get_W(pxx, N, M, W);

    if (DEBUG) Rprintf("Compute information matrix\n");
    get_Info2(tXXZ, N, M, Nparm, W, W_index, Info);

    if (DEBUG) Rprintf("Compute covariance matrix\n");
    rc = cov_inv(Info, Znc, Inv);
    if (rc) {
      Rprintf("ERROR computing inverse of information matrix\n");
      error("ERROR");
    }

    if (DEBUG) Rprintf("Compute W_y\n");
    get_Wy2(Y, lxx, pxx, W, W_index, N, M, w_y);

    if (DEBUG) Rprintf("Compute delta\n");
    get_delta(Inv, Nparm, tXXZ, N, M, w_y, ret_delta);

    if (DEBUG > 1) print_dVec(ret_delta, Nparm, "delta");

    /* Check for non-finite values */
    if (!all_finite(ret_delta, Nparm)) {
      Rprintf("ERROR: algorithm not converging, parameters have non-finite values\n");
      error("ERROR");
    }

    rerror = checkStop(ret_delta, delta0, Nparm);
    if (DEBUG) Rprintf("Check stopping criteria: relative error = %g\n", rerror);
    if (rerror <= tol) {
      conv = 1;
      break;
    }

    /* Update delta0 */
    if (DEBUG) Rprintf("Update parameters\n");
    for (i=0; i<Nparm; i++) delta0[i] = ret_delta[i];
  }

  if (DEBUG) Rprintf("Free memory\n");
  free(lxx);
  matrix_free((void **)Info, Nparm);
  free(w_y);

  if (!conv) {
    Rprintf("ERROR: algorithm did not converge\n");
    error("ERROR");
  }

} /* END: mvpoly2 */

/* Function to compute t(XXZ_test) */
static void tXXZvec(Xvec, N, Z, M, Znc, out)
double *Xvec, **Z, **out; /* out must have dim Znc X NM */
int N, M, Znc;
{
  int i, j, k;
  double **pro, *pco, *px, zi;

  for (i=0, pro=out; i<Znc; i++, pro++) {
    pco = *pro;
    for (j=0; j<M; j++) {
      zi = Z[j][i];
      for (k=0, px=Xvec; k<N; k++, px++) {
        *pco++ = *px * zi;
      }
    }
  }

} /* END: tXXZvec */

/* Function to compute the score tXXZ*(Y - mu) */
static void get_score(tXXZ, Y, Mu, N, M, Znc, out)
double **tXXZ, *Y, *Mu, *out;
int N, M, Znc;
{
  int i, j, NM;
  double *po, *py, *pmu, **prx, *pcx, sum;

  NM = N*M;
  for (i=0, po=out, prx=tXXZ; i<Znc; i++, po++, prx++) {
    sum = 0.0;
    for (j=0, py=Y, pmu=Mu, pcx=*prx; j<NM; j++, py++, pmu++, pcx++) {
      sum += *pcx * (*py - *pmu); 
    }
    *po = sum;
  }

} /* END: get_score */

/* Function to compute information matrix t(xxz1)%*%W%*%xxz2 */
static void get_quadFormW(tXXZ1, tXXZ2, W, W_index, N, M, nc1, nc2, out)
double **tXXZ1, **tXXZ2, **out, **W;
int N, M, nc1, nc2, **W_index;
{
  int i, j, row, col, NM, **pr, *pc, cnt;
  double sum, **prW, *pcW, sum2, **pr1, *pc1, *pc2;

  NM = N*M;
  pr = W_index;

  for (row=0, pr1=tXXZ1; row<nc1; row++, pr1++) {
    for (col=0; col<nc2; col++) {
      /* Compute tXXZ1[row, ]*W*XXZ2[, col] = tXXZ1[row, ]*W*tXXZ2[col, ] */
      sum2 = 0.0;
      cnt  = 1;
      pc2  = tXXZ2[col];
      for (i=0, prW=W, pc1=*pr1; i<NM; i++, prW++, pc1++) {
        sum = 0.0;
        for (j=0, pc=*pr, pcW=*prW; j<M; j++, pc++, pcW++) {
          sum += *pcW * pc2[*pc];
        } 
        sum2 += *pc1 * sum;
        cnt++;
        if (cnt > N) {
          cnt = 1;
          pr  = W_index;
        } else {
          pr++;
        }
      }
      out[row][col] = sum2;
    }
  }

} /* END: get_quadFormW */

/* Function to compute a quadratic form t(X)*A*X */
static void quadForm(X, A, nrx, ncx, out)
double **X, **A, **out;
int nrx, ncx;
{
  int row, col, i, j;
  double sum, sum2;

  for (row=0; row<ncx; row++) {
    for (col=row; col<ncx; col++) {
      sum2 = 0.0;
      for (i=0; i<nrx; i++) {
        sum = 0.0;
        for (j=0; j<nrx; j++) {
          sum += A[i][j]*X[j][col];
        }
        sum2 += X[i][row]*sum;
      }
      out[row][col] = sum2;
      if (row != col) out[col][row] = sum2; 
    }
  }

} /* END: quadForm */

/* Function to compute a quadratic form with a vector */
static double quadFormVec(x, A, n)
double *x, **A;
int n;
{
  int i, j;
  double **pr, *pc, *px, *px2, sum, sum2;

  sum2 = 0.0;
  for (i=0, px=x, pr=A; i<n; i++, px++, pr++) {
    sum = 0.0;
    for (j=0, px2=x, pc=*pr; j<n; j++, px2++, pc++) sum += *pc * *px2;
    sum2 += *px * sum; 
  }
  
  return(sum2);

} /* END: quadFormVec */

/* Function to compute information matrix in score_test function */
static void getInfoS(tXXZ_test, W, W_index, tXXZ, Inv, N, M, Nparm, Ncatp1, out)
double **tXXZ_test, **W, **tXXZ, **Inv,  **out;
int **W_index, N, M, Nparm, Ncatp1;
{
  int row, col;
  double **Q, **Q0, **tmp, val;

  Q   = dMat_alloc(Nparm, Ncatp1, 0, 0.0);
  Q0  = dMat_alloc(Ncatp1, Ncatp1, 0, 0.0);
  tmp = dMat_alloc(Ncatp1, Ncatp1, 0, 0.0);

  /* Compute Q = tXXZ*W*XXZ_test */ 
  get_quadFormW(tXXZ, tXXZ_test, W, W_index, N, M, Nparm, Ncatp1, Q);

  /* Compute Q0 = tXXZ_test*W*XXZ_test. This can be made more efficient. */ 
  get_quadFormW(tXXZ_test, tXXZ_test, W, W_index, N, M, Ncatp1, Ncatp1, Q0);

  quadForm(Q, Inv, Nparm, Ncatp1, tmp);

  for (row=0; row<Ncatp1; row++) {
    for (col=row; col<Ncatp1; col++) {
      val = Q0[row][col] - tmp[row][col];
      out[row][col] = val;
      if (row != col) out[col][row] = val;  
    }
  }

  matrix_free((void **)Q, Nparm);
  matrix_free((void **)Q0, Ncatp1);
  matrix_free((void **)tmp, Ncatp1);

} /* END: getInfoS */ 

/* Function to compute the p-value */
static double chiTestP(score, inv, n)
double *score, **inv;
int n;
{
  double test, p;

  test = quadFormVec(score, inv, n);
  p    = pchisq(test, (double) n, 0, 0);
  
  return(p);

} /* END: chiTestP */

/* Function score_test */
void score_test(deltai, pNparm, Y, Xvec, Xtest, Zvec, pN, pM, pNcat, pNcov, pNiter, ptol, 
            pDEBUG, ret_rc, ret_delta, ret_p)
double *deltai, *Y, *Xvec, *Xtest, *Zvec, *ptol, *ret_delta, *ret_p;
int *pNparm, *pN, *pM, *pNcat, *pNcov, *pNiter, *ret_rc, *pDEBUG;
{
  int Niter, M, N, Ncat, Ncatp1, Ncov0, Ncov, Znr, Znc, NM;
  int Nparm, DEBUG, **W_index;
  double tol, **X, **Z_design, **Z, **W, *score;
  double **Inv, **InfoS, *pxx, **tXXZ_test, **tXXZ, pval;

  *ret_rc = 1;
  DEBUG   = *pDEBUG;
  Nparm   = *pNparm;
  Niter   = *pNiter;
  N       = *pN;
  M       = *pM;
  Ncat    = *pNcat;
  Ncatp1  = Ncat + 1;
  tol     = *ptol;
  Ncov0   = *pNcov;
  Ncov    = Ncov0 + 1;  /* Allow for intercept */
  Znr     = M*Ncov;
  Znc     = (Ncov-1)*Ncatp1+M;
  NM      = N*M;
 
  if (Nparm != Znc) error("Nparm != Znc\n");

  /* Allocate memory for matrix of covariates and Z map matrix */
  if (DEBUG) Rprintf("Allocate memory\n");
  Z_design  = dMat_alloc(M, Ncatp1, 0, 0.0);
  Z         = dMat_alloc(Znr, Znc, 1, 0.0); /* Initialize to 0 */
  W         = dMat_alloc(NM, M, 0, 0.0);
  W_index   = iMat_alloc(N, M, 0, 0);
  tXXZ_test = dMat_alloc(Ncatp1, NM, 0, 0.0);
  tXXZ      = dMat_alloc(Nparm, NM, 0, 0.0);
  X         = dMat_alloc(N, Ncov, 0, 0.0);
  pxx       = dVec_alloc(NM, 0, 0.0);
  Inv       = dMat_alloc(Nparm, Nparm, 0, 0.0);
  InfoS     = dMat_alloc(Ncatp1, Ncatp1, 0, 0.0);
  score     = dVec_alloc(Ncatp1, 0, 0.0);
  printf("%i\n",DEBUG);

  if (DEBUG) Rprintf("Copy data\n");
  printf("%i\n",DEBUG);
  fillMat(Xvec, N, Ncov0, 1, X);
  /*print_dMat(X,N,Ncov,X); */
  fillMat(Zvec, M, Ncatp1, 0, Z_design);
  /*printf("N\tM\tNcat\tNcatp1\tNcov0\tNcov\tZnr\t\Znc\n");
  printf("%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",N,M,Ncat,Ncatp1,Ncov0,Ncov,Znr,Znc);*/
    /*printf("%i\n",DEBUG);*/
  
    if (DEBUG) Rprintf("Get the matrix Z\n");
    get_Z(Z_design, M, Ncatp1, Ncov, Z);
  /*print_dMat(Z,Znr,Znc,Z);*/ 
  if (DEBUG) Rprintf("Get the matrix W_index\n");
  get_W_index(N, M, W_index);
  if (DEBUG) Rprintf("Get the matrix t(XXZ)\n");
  get_tXXZ(X, N, Ncov, M, Ncatp1, Z, tXXZ);
  if (DEBUG) Rprintf("Call mvpoly2\n");
  mvpoly2(deltai, Nparm, Y, X, Z, N, M, Ncat, Ncov, Niter, tol, 
            ret_delta, W, W_index, Inv, pxx, tXXZ, DEBUG);
  matrix_free((void **)X, N);
  matrix_free((void **)Z, M);
  
  if (DEBUG) Rprintf("Call tXXZvec\n");
  tXXZvec(Xtest, N, Z_design, M, Ncatp1, tXXZ_test);

  if (DEBUG) Rprintf("Compute the score\n");
  get_score(tXXZ_test, Y, pxx, N, M, Ncatp1, score);
  free(pxx);

  if (DEBUG) Rprintf("Compute the information matrix\n");
  getInfoS(tXXZ_test, W, W_index, tXXZ, Inv, N, M, Nparm, Ncatp1, InfoS);
  matrix_free((void **) W, NM);
  matrix_free((void **) W_index, N);
  matrix_free((void **)tXXZ_test, Ncatp1);
  matrix_free((void **)tXXZ, Nparm);
  matrix_free((void **)Z_design, M);
  matrix_free((void **)Inv, Nparm);

  /* Compute inverse */
  if (DEBUG) Rprintf("Compute the inverse of the information matrix\n");
  Inv = dMat_alloc(Ncatp1, Ncatp1, 0, 0.0);
  cov_inv(InfoS, Ncatp1, Inv);

  if (DEBUG) Rprintf("Compute final p-value\n");
  pval = chiTestP(score, Inv, Ncatp1);
  *ret_p = pval;

  if (DEBUG) Rprintf("Free memory\n");
  matrix_free((void **)InfoS, Ncatp1);
  matrix_free((void **)Inv, Ncatp1);
  free(score);

  *ret_rc = 0;
  return;

} /* END: score_test */




