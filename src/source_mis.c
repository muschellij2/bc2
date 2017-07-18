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

}

  /* Get the column for diag(P) (X) z_design, loop over rows of z_design */
  /*for(i=0; i<M; i++) {
    vec = Z_design[i];
    for (j=0; j<Ncov; j++) {
      col = j*Ncatp1;
      for (k=0; k<Ncatp1; k++) {
        out[row][col] = vec[k];
        col++;
      }
      row++;
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





/* Function to compute information matrix t(xxz)%*%W%*%xxz */
static void get_Info(tXXZ, N, M, nparm, Pxx, out)
double **tXXZ, *Pxx, **out;
int N, M, nparm;
{
  int row, col;
  double val, *prow, *pcol;

  for (row=0; row<nparm; row++) {
    prow = tXXZ[row];
    for (col=row; col<nparm; col++) {
      /* Compute XXZ[, row]*W*XXZ[, col] */
      val = v1Wv2(Pxx, N, M, prow, tXXZ[col]);
      out[row][col] = val;
      if (row != col) out[col][row] = val; 
    }
  }

} /* END: get_Info */

 /* fill the info matrix to the result*/
static void fill_Info(Info,ret_info,Nparm)
double **Info,*ret_info;
int Nparm;
{
  int i,j;
for(j=0;j<Nparm;j++){
  for(i=0;i< Nparm;i++){
    ret_info[i*Nparm+j] = Info[i][j];
  }
}


}

/*end fill_Info*/





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




  /* Zero out the diagonal and above */


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
 /* END: symPosMatInv */

 /* END: cov_inv */
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


/* Function Mvpoly */
void mis_infor(pNparm, Y, Xvec, Zvec, pN, pM, pNcat, pNcov, pNiter, ptol, 
            pDEBUG, ret_rc,ret_info)
double *Y, *Xvec, *Zvec, *ptol,*ret_info;
int *pNparm, *pN, *pM, *pNcat, *pNcov, *pNiter, *ret_rc, *pDEBUG;
{
  int i, Niter, M, N, Ncat, Ncatp1, Ncov0, Ncov, Znr, Znc, NM, rc;
  int Nparm, DEBUG;
  double tol, **X, **Z_design, **Z;
  double  **Info, **tXXZ;

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
 /*printf("Ncat\tNcatp1\tNcov0\tNcov\tZnr\Znc\n");
 printf("%i\t%i\t%i\t%i\t%i\t%i\n",Ncat,Ncatp1,Ncov0,Ncov,Znr,Znc);*/
  

  /* Allocate memory for matrix of covariates and Z map matrix */

  tXXZ     = dMat_alloc(Nparm, NM, 0, 0.0);

  X        = dMat_alloc(N, Ncov, 0, 0.0);
  Z_design = dMat_alloc(M, Ncatp1, 0, 0.0);
  Z        = dMat_alloc(Znr, Znc, 1, 0.0); /* Initialize to 0 */
  Info     = dMat_alloc(Nparm, Nparm, 0, 0.0);
 
 
  /* Copy initial estimates to delta0 */

  fillMat(Xvec, N, Ncov0, 1, X);
  fillMat(Zvec, M, Ncatp1, 0, Z_design);
  get_Z(Z_design, M, Ncatp1, Ncov, Z);
  /*print_dMat(Z,Ncatp1*M,(Ncov-1)*Ncatp1+M,Z);*/

  /* Get the matrix t(XXZ) */
 
  get_tXXZ(X, N, Ncov, M, Ncatp1, Z, tXXZ);
  /*print_dMat(tXXZ,(Ncov-1)*Ncatp1+M,M*N,tXXZ);*/

  

    /*get_Info(X, N, Ncov, M, Z, Znr, Znc, pxx, Info);*/
    get_Info(tXXZ, N, M, Nparm, Y, Info);

  matrix_free((void **)tXXZ, Nparm);
  matrix_free((void **)X, N);
  matrix_free((void **)Z_design, M);
  matrix_free((void **)Z, M);
  fill_Info(Info,ret_info,Nparm);
  matrix_free((void **)Info, Znc);


  return;

} /* END: Mvpoly */

