#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#include "protos.h"

#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

int qsplitC(FLOAT *a, int *ind, int n, int ncut)
{
/*----------------------------------------------------------------------
|     does a quick-sort split of a complex real array.
|     on input a[0 : (n-1)] is a real array
|     on output is permuted such that its elements satisfy:
|
|     abs(a[i]) >= abs(a[ncut-1]) for i < ncut-1 and
|     abs(a[i]) <= abs(a[ncut-1]) for i > ncut-1
|
|     ind[0 : (n-1)] is an integer array permuted in the same way as a.
|---------------------------------------------------------------------*/
   FLOAT tmp;
   double abskey;
   int j, itmp, first, mid, last;
   first = 0;
   last = n-1;
   if (ncut<first || ncut>last) return 0;
/* outer loop -- while mid != ncut */
label1:
   mid = first;
   abskey = ABS_VALUE(a[mid]);
  for (j=first+1; j<=last; j++) {
     if (ABS_VALUE(a[j]) > abskey) {
	 tmp = a[++mid];
	 itmp = ind[mid];
	 a[mid] = a[j];
	 ind[mid] = ind[j];
	 a[j]  = tmp;
	 ind[j] = itmp;
      }
   }
/*-------------------- interchange */
   tmp = a[mid];
   a[mid] = a[first];
   a[first]  = tmp;
   itmp = ind[mid];
   ind[mid] = ind[first];
   ind[first] = itmp;
/*-------------------- test for while loop */
   if (mid == ncut) return 0;
   if (mid > ncut) 
      last = mid-1;
   else
      first = mid+1;
   goto label1;
}
/*--------------- end of zqsplitC ----------------------------------------*/

int SparTran(csptr amat, csptr bmat, int job, int flag)
{
/*----------------------------------------------------------------------
| Finds the transpose of a matrix stored in SparRow format.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SparRow format.
|
| job    = integer to indicate whether to fill the values (job.eq.1)
|          of the matrix (bmat) or only the pattern.
|
| flag   = integer to indicate whether the matrix has been filled
|          0 - no filled
|          1 - filled
|
| on return:
| ----------
| (bmat) = the transpose of (mata) stored in SparRow format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
  int i, j, *ind, pos, size=amat->n, *aja;
  FLOAT *ama=NULL;
  ind = (int *) Malloc(size*sizeof(int), "SparTran:1" );
  for (i=0; i<size; i++)
    ind[i] = 0;
  if(!flag) {
/*--------------------  compute lengths  */
    for (i=0; i<size; i++) {
      aja = amat->pj[i];
      for (j=0; j<amat->nnzrow[i]; j++)
	ind[aja[j]]++;
    }
/*--------------------  allocate space  */
    for (i=0; i<size; i++) {
      bmat->pj[i] = (int *) Malloc(ind[i]*sizeof(int), "SparTran:2" );
      bmat->nnzrow[i] = ind[i];
      if (job == 1) {
	bmat->pa[i] = (FLOAT *) Malloc(ind[i]*sizeof(FLOAT), "SparTran:3" );
      }
      ind[i] = 0;
    }
  }
/*--------------------  now do the actual copying  */
  for (i=0; i<size; i++) {
    aja = amat->pj[i];
    if (job == 1)
      ama = amat->pa[i];
    for (j=0; j<amat->nnzrow[i]; j++) {
      pos = aja[j];
      bmat->pj[pos][ind[pos]] = i;
      if (job == 1)
	bmat->pa[pos][ind[pos]] = ama[j];
      ind[pos]++;
    }
  }
  free(ind);
  return 0;
}
/*-------------- end of SparTran ---------------------------------------
|---------------------------------------------------------------------*/

void swapj(int v[], int i, int j)
{
  int temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}
void swapm(FLOAT v[], int i, int j) 
{
  FLOAT temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}

void dswapm(double v[], int i, int j)
{
   double temp;
   
   temp = v[i];
   v[i] = v[j];
   v[j] = temp;
}

int roscalC(csptr mata, double *diag, int nrm)
{
/*---------------------------------------------------------------------
|
| This routine scales each row of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SparRow form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> row j is a zero row
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, k;
   FLOAT *kr; 
   double scal;

   for (i=0; i<mata->n; i++) {
      scal = 0.0;
      kr = mata->pa[i];
      if (nrm == 0) {
	 for (k=0; k<mata->nnzrow[i]; k++)
	    if (ABS_VALUE(kr[k]) > fabs(scal)) scal = ABS_VALUE(kr[k]);
      }
      else if (nrm == 1) {
         for (k=0; k<mata->nnzrow[i]; k++)
            scal += ABS_VALUE(kr[k]);
      }
      else {  /* nrm = 2 */
         for (k=0; k<mata->nnzrow[i]; k++)
            scal += ABS_VALUE(kr[k]*kr[k]);
      }
      if (nrm == 2) scal = sqrt(scal);
      if(ABS_VALUE(scal-DBL_EPSILON) <= DBL_EPSILON*ABS_VALUE(scal)){
	scal = 1.0; 
	/* YS. return i+1; */
      }
      else 
	scal = 1.0 / scal;
      diag[i] = scal;
      for (k=0; k<mata->nnzrow[i]; k++)
	 kr[k] = kr[k] * scal;
   }
   return 0;
}
/*---------------end of roscalC-----------------------------------------
----------------------------------------------------------------------*/
int coscalC(csptr mata, double *diag, int nrm)
{
/*---------------------------------------------------------------------
|
| This routine scales each column of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SparRow form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> column j is a zero column
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, j, k;
   FLOAT *kr;
   int *ki;

   for (i=0; i<mata->n; i++)
      diag[i] = 0.0;
/*---------------------------------------
|   compute the norm of each column
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      kr = mata->pa[i];
      ki = mata->pj[i];
      if (nrm == 0) {
	 for (k=0; k<mata->nnzrow[i]; k++) {
	    j = ki[k];
	    if (ABS_VALUE(kr[k]) > diag[j]) diag[j] = ABS_VALUE(kr[k]);
	 }
      }
      else if (nrm == 1) {
         for (k=0; k<mata->nnzrow[i]; k++)
            diag[ki[k]] += ABS_VALUE(kr[k]);
      }
      else {  /*  nrm = 2 */
         for (k=0; k<mata->nnzrow[i]; k++)
            diag[ki[k]] += ABS_VALUE(kr[k]*kr[k]);
      }
   }
   if (nrm == 2) {
      for (i=0; i<mata->n; i++)
	 diag[i] = sqrt(diag[i]);
   }
/*---------------------------------------
|   invert
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      if(ABS_VALUE(diag[i]-DBL_EPSILON) <= DBL_EPSILON*ABS_VALUE(diag[i]))
	/* return i+1;*/
	diag[i] = 1.0; 
      else 
	 diag[i] = 1.0 / diag[i];
   }
/*---------------------------------------
|   C = A * D
|--------------------------------------*/
   for (i=0; i<mata->n; i++) {
      kr = mata->pa[i];
      ki = mata->pj[i];
      for (k=0; k<mata->nnzrow[i]; k++)
	 kr[k] = kr[k] * diag[ki[k]];
   }
   return 0;
}
/*---------------end of coscalC-----------------------------------------
----------------------------------------------------------------------*/
void dscale(int n, double *dd, FLOAT *x, FLOAT * y)
{ 
/* Computes  y == DD * x                               */
/* scales the vector x by the diagonal dd - output in y */
  int k;

  for (k=0; k<n; k++) 
    y[k] = dd[k]*x[k];
}

  
void qsortC(int *ja, FLOAT *ma, int left, int right, int abval){
/*----------------------------------------------------------------------
|
| qqsort: sort ma[left]...ma[right] into decreasing order
| from Kernighan & Ritchie
|
| ja holds the column indices
| abval = 1: consider absolute values
|         0: values
|
|---------------------------------------------------------------------*/
  int i, last;

  if (left >= right)  return;
  if (abval) {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (ABS_VALUE(ma[i]) > ABS_VALUE(ma[left])) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsortC(ja, ma, left, last-1, abval);
    qsortC(ja, ma, last+1, right, abval);
  }
  else {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (ABS_VALUE(ma[i]) > ABS_VALUE(ma[left])) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsortC(ja, ma, left, last-1, abval);
    qsortC(ja, ma, last+1, right, abval);
  }
}

void printmat(FILE *ft, csptr A, int i0, int i1){
/*-------------------------------------------------------------+
| to dump rows i0 to i1 of matrix for debugging purposes       |
|--------------------------------------------------------------*/
  int i, k, nzi;
  int *row;
  FLOAT *rowm;
  for (i=i0; i<i1; i++)    {
    nzi = A->nnzrow[i];
    row = A->pj[i];
    rowm = A->pa[i];
#if defined(DBL_CMPLX)    
    for (k=0; k< nzi; k++){
      fprintf(ft," row %d  a  (%e, %e) ja %d \n", i+1, creal(rowm[k]), cimag(rowm[k]), row[k]+1);
    }
#else
    for (k=0; k< nzi; k++){
      fprintf(ft," row %d  a  %e ja %d \n", i+1, rowm[k], row[k]+1);
    }
#endif
  }
}

void qsortR2I(double *wa, int *cor1, int *cor2, int left, int right){
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into decreasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
  int i, last;
/*
  void swapj(int *, int, int);
  void dswapm(double *, int, int);
*/  
  if (left >= right)  return;
  
  dswapm(wa, left, (left+right)/2);
  swapj(cor1, left, (left+right)/2);
  swapj(cor2, left, (left+right)/2);
  last = left;
  for (i=left+1; i<=right; i++) {
    if (wa[i] > wa[left]) {
      dswapm(wa, ++last, i);
      swapj(cor1, last, i);
      swapj(cor2, last, i);
    }
  }
  dswapm(wa, left, last);
  swapj(cor1, left, last);
  swapj(cor2, left, last);
  qsortR2I(wa, cor1, cor2, left, last-1);
  qsortR2I(wa, cor1, cor2, last+1, right);
}


void qsort2C(int *ja, FLOAT *ma, int left, int right, int abval){
/*----------------------------------------------------------------------
|
| qqsort: sort ma[left]...ma[right] into increasing order
| from Kernighan & Ritchie
|
| ja holds the column indices
| abval = 1: consider absolute values
|         0: values
|
|---------------------------------------------------------------------*/
  int i, last;
  if (left >= right)  return;
  if (abval) {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (ABS_VALUE(ma[i]) < ABS_VALUE(ma[left])) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsort2C(ja, ma, left, last-1, abval);
    qsort2C(ja, ma, last+1, right, abval);
  }
  
  else {
    swapj(ja, left, (left+right)/2);
    swapm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
      if (ABS_VALUE(ma[i]) < ABS_VALUE(ma[left])) {
	swapj(ja, ++last, i);
	swapm(ma, last, i);
      }
    }
    swapj(ja, left, last);
    swapm(ma, left, last);
    qsort2C(ja, ma, left, last-1, abval);
    qsort2C(ja, ma, last+1, right, abval);
  }
}

void qqsort(int *ja, FLOAT *ma, int left, int right){
/*----------------------------------------------------------------------
|
| qqsort: sort ja[left]...ja[right] into increasing order
| from Kernighan & Ritchie
|
| ma holds the real values
|
|---------------------------------------------------------------------*/
  int i, last;
  if (left >= right)  return;
  swapj(ja, left, (left+right)/2);
  swapm(ma, left, (left+right)/2);
  last = left;
  for (i=left+1; i<=right; i++) {
    if (ja[i] < ja[left]) {
      swapj(ja, ++last, i);
      swapm(ma, last, i);
    }
  }
  swapj(ja, left, last);
  swapm(ma, left, last);
  qqsort(ja, ma, left, last-1);
  qqsort(ja, ma, last+1, right);
}

void hilosort(csptr mat, int abval, int hilo){
/*----------------------------------------------------------------------
|
| This routine sorts the entries in each row of a matrix from hi to low.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SparRow format.
|
| abval =   1: use absolute values of entries
|           0: use values
|
| hilo  =   1: sort in decreasing order
|           0: sort in increasing order
|
|
| on return:
| ----------
| (mat) = (mat) where each row is sorted.
|
|---------------------------------------------------------------------*/
  int j, n=mat->n, *nnz=mat->nnzrow;
/*
  void qsortC(int *, FLOAT *, int, int, int);
  void qsort2C(int *, FLOAT *, int, int, int);
*/  
  if (hilo)
    for (j=0; j<n; j++)
      qsortC(mat->pj[j], mat->pa[j], 0, nnz[j]-1, abval);
  
  else
    for (j=0; j<n; j++)
      qsort2C(mat->pj[j], mat->pa[j], 0, nnz[j]-1, abval);
  
  return;
}
/*------- end of hilosort ----------------------------------------------
|---------------------------------------------------------------------*/

void qsort3i(int *wa, int *cor1, int *cor2, int left, int right)
/*----------------------------------------------------------------------
|
| qqsort: sort wa[left]...wa[right] into increasing order
| from Kernighan & Ritchie
|
|---------------------------------------------------------------------*/
{
   int i, last;
/*
   void swapj(int *, int, int);
*/
   if (left >= right)  return;

   swapj(wa, left, (left+right)/2);
   swapj(cor1, left, (left+right)/2);
   swapj(cor2, left, (left+right)/2);
   last = left;
   for (i=left+1; i<=right; i++) {
      if (wa[i] < wa[left]) {
	 swapj(wa, ++last, i);
	 swapj(cor1, last, i);
	 swapj(cor2, last, i);
      }
   }
   swapj(wa, left, last);
   swapj(cor1, left, last);
   swapj(cor2, left, last);
   qsort3i(wa, cor1, cor2, left, last-1);
   qsort3i(wa, cor1, cor2, last+1, right);
}

