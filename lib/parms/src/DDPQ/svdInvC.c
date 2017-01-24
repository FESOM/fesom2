#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "protos.h"

#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

#define TOL 1.e-17
#define max_abs(a,b) ((ABS_VALUE(a)>ABS_VALUE(b))?(a):(b))

int invGauss(int nn, FLOAT *A) 
{
  /* *-------------------- inversion by svd
     This calls lapack routines for inverting a dense matrix.
     dgetrf and dgetri

     ON ENTRY
     ** A = square matrix of size n x n -- dimensioned with
     ** leading dimension =  n

     ON RETURN A contains the inverse  of the input matrix.
  */
  int lWk, info;

  FLOAT *Wk;
  int *ipiv;

  lWk = 10*nn;

  /*-------------------- trivial case nn = 1                     */
  if (nn == 1) {
    if(ABS_VALUE(A[0]-DBL_EPSILON) <= DBL_EPSILON*ABS_VALUE(A[0]))
      return 1;
    else {
      A[0] = 1.0 / A[0];
      return 0;
    }
  }
  /*-------------------- general case                              */

  Wk  = (FLOAT *) malloc(lWk*sizeof(FLOAT));
  ipiv = (int *) malloc(nn*sizeof(int));
  if (Wk == NULL || ipiv == NULL)
    return -1;
  /*-------------------- get LU factorization with pivoting         */
  GGETRF (nn, nn,A, nn, ipiv, info);  

  if (info !=0 ) return info;
  /*-------------------- compute inverse ---- */     
  GGETRI (nn, A, nn, ipiv, Wk, lWk, info);

  free(Wk);
  free(ipiv);
  return info;
}

int invSVD(int nn, FLOAT *A) 
{
  /* *-------------------- inversion by svd
     This calls lapack routine dgesvd --
     ON ENTRY
     ** A = square matrix of size n x n -- dimensioned with
     ** leading dimension =  n
     ON RETURN A contains the truncated SVD inverse of input matrix.
     ** tolerance set for truncation is TOL and can be changed in
     ** above define statement
     **--------------------
     */
  int lWk, i, info;

  double *S, *Rwk;
  FLOAT *U, *VT, *Wk;
  REAL nrm;
  FLOAT tmp, one=1.0, zero=0.0;

  double tol=TOL;

  lWk = 5*nn;

  U  = (FLOAT *) malloc(nn*nn*sizeof(FLOAT));
  VT = (FLOAT *) malloc(nn*nn*sizeof(FLOAT));
  S  = (double *) malloc(nn*sizeof(double));
  Rwk  = (double *) malloc(5*nn*sizeof(double));  
  Wk  = (FLOAT *) malloc(lWk*sizeof(FLOAT));

  if (U == NULL || VT == NULL || S == NULL || Wk == NULL)
    return -1;
  /*-------------------- trivial case nn = 1                     */
  if (nn == 1) {
    if(ABS_VALUE(A[0]-DBL_EPSILON) <= DBL_EPSILON*ABS_VALUE(A[0]))
      return 1;
    else {
      A[0] = one / A[0];
      return 0;
    }
  }
  /*-------------------- general case                              */
#if defined(DBL_CMPLX)
  zgesvd_ ("A","A", &nn, &nn, A, &nn, S, U, &nn, VT, &nn, Wk, &lWk,
	  Rwk, &info) ;
#else
  dgesvd_ ("A","A", &nn, &nn, A, &nn, S, U, &nn, VT, &nn, Wk, &lWk,
	  &info) ;
#endif	  
	  
  if(ABS_VALUE(S[0]-DBL_EPSILON) <= DBL_EPSILON*ABS_VALUE(S[0]))
    return 1;

  nrm = S[0]*tol;
  /*-------------------- compute S\inv * VT                        */
  for (i=0; i<nn; i++) {
    tmp = one / max_abs(S[i],nrm) ;
    GSCAL(nn, tmp,&VT[i],nn);
  }
  /*-------------------- do [V^T S\inv ] * U^T                     */
  GGEMM("t","t",nn,nn,nn, one, VT, nn, U, nn, zero, A, nn);
  /*-------------------- Done -------------------------------------*/
  free(U);
  free(VT);
  free(S);
  free(Wk);
  free(Rwk);
  return 0;
}

