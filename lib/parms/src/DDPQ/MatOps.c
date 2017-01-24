#include <stdio.h>
#include <stdlib.h>
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#include "protos.h"

void matvec( csptr mata, FLOAT *x, FLOAT *y )  
{
  /*---------------------------------------------------------------------
    | This function does the matrix vector product y = A x.
    |----------------------------------------------------------------------
    | on entry:
    | mata  = the matrix (in SparRow form)
    | x     = a vector
    |
    | on return
    | y     = the product A * x
    |--------------------------------------------------------------------*/
  /*   local variables    */
  int i, k, *ki;
  FLOAT *kr;

  for (i=0; i<mata->n; i++) {
    y[i] = 0.0;
    kr = mata->pa[i];
    ki = mata->pj[i];
    for (k=0; k<mata->nnzrow[i]; k++)
      y[i] += kr[k] * x[ki[k]];
  }
  return;
}

void invsp(int start, ilutptr ilusch, FLOAT *y, FLOAT *x)
{
  /*---------------------------------------------------------------------
    |
    | This routine does the backward solve U x = y, where U is upper or
    | bottom part of local upper triangular matrix
    |
    | Can be done in place.
    |
    | Zhongze Li, Aug. 17th, 2001
    |
    |----------------------------------------------------------------------
    | on entry:
    | start = the index of the first component
    | n     = one ore than the index owned by the processor 
    | y     = a vector
    | ilusch  = the LU matrix as provided from the ILU routines.
    |
    | on return
    | x     = the product U^{-1} * y
    |
    |---------------------------------------------------------------------*/
  /*   local variables    */
  int i, k, *ki, n;
  FLOAT *kr, t;

  n = ilusch->L->n;

  for (i=start; i < n; i++) {
    t = y[i-start];
    if ( ilusch->L->nnzrow[i] > 0 ) {
      kr = ilusch->L->pa[i];
      ki = ilusch->L->pj[i];
      for (k=0; k<ilusch->L->nnzrow[i]; k++)
	if(ki[k] >= start && ki[k] < n) {
	  t -= kr[k]*y[ki[k]-start];
	}
    }
    x[i-start] = t;
  }

  for (i=n-1; i>= start; i--) {
    kr = ilusch->U->pa[i];
    ki = ilusch->U->pj[i];
    t  = x[i-start];
    for (k=1; k<ilusch->U->nnzrow[i]; k++)
      t -= kr[k] * x[ki[k]-start];
    x[i-start] = t*kr[0];
  }
  return;
}

void arms_Lsol(csptr mata, FLOAT *b, FLOAT *x)
{
  /*---------------------------------------------------------------------
    | This function does the forward solve L x = b.
    | Can be done in place.
    |----------------------------------------------------------------------
    | on entry:
    | mata  = the matrix (in SparRow form)
    | b     = a vector
    |
    | on return
    | x     = the solution of L x = b 
    |--------------------------------------------------------------------*/
  /*   local variables    */
  int i, k;
  FLOAT *kr;
  int *ki;
  for (i=0; i<mata->n; i++) {
    x[i] = b[i];
    if ( mata->nnzrow[i] > 0 ) {
      kr = mata->pa[i];
      ki = mata->pj[i];
      for (k=0; k<mata->nnzrow[i]; k++)
	x[i] -= kr[k]*x[ki[k]];
    }
  }
  return;
}
/*---------------end of arms_Lsol-----------------------------------------
----------------------------------------------------------------------*/
void arms_Usol(csptr mata, FLOAT *b, FLOAT *x)
{
  /*---------------------------------------------------------------------
    | This function does the backward solve U x = b.
    | Can be done in place.
    |----------------------------------------------------------------------
    | on entry:
    | mata  = the matrix (in SparRow form)
    | b    = a vector
    |
    | on return
    | x     = the solution of U * x = b 
    |
    |---------------------------------------------------------------------*/
  /*   local variables    */
  int i, k, *ki;
  FLOAT *kr;
  for (i=mata->n-1; i>=0; i--) {
    kr = mata->pa[i];
    ki = mata->pj[i];
    x[i] = b[i] ;
    for (k=1; k<mata->nnzrow[i]; k++)
      x[i] -= kr[k] * x[ki[k]];
    x[i] *= kr[0];
  }
  return;
}
/*----------------end of arms_Usol----------------------------------------
----------------------------------------------------------------------*/

int descend(p4ptr levmat, FLOAT *x, FLOAT *wk)
{
/*---------------------------------------------------------------------
| This function does the (block) forward elimination in ARMS
|                       new       old
|     |            |  |     |    |    |
|     | L        0 |  | wx1 |    | x1 |
|     |            |  |     | =  |    | 
|     | EU^{-1}  I |  | wx2 |    | x2 |
|     |            |  |     |    |    |
| x used and not touched -- or can be the same as wk.
|--------------------------------------------------------------------*/
/*  local variables   */
  int j, len=levmat->n, lenB=levmat->nB, *iperm=levmat->rperm; 
  FLOAT *work = levmat->wk; 
/*------------------------------------------------------
|   apply permutation P to rhs 
|-----------------------------------------------------*/
  for (j=0; j<len; j++)
    work[iperm[j]] = x[j] ;
  arms_Lsol(levmat->L, work, wk);       /* sol:   L x = x                 */
  arms_Usol(levmat->U, wk, work);        /* sol:   U work(2) = work         */
/*-------------------- compute x[lenb:.] = x [lenb:.] - E * work(1) */
  matvecz (levmat->E, work, &work[lenB], &wk[lenB]) ; 
  return 0;
}
/*----end-of-descend---------------------------------------------------
|----------------------------------------------------------------------
|--------------------------------------------------------------------*/
int ascend (p4ptr levmat, FLOAT *x, FLOAT *wk) 
{
  /*---------------------------------------------------------------------
    | This function does the (block) backward substitution: 
    |
    |     |            |  |     |    |    |
    |     | U  L^{-1}F |  | wk1 |    | x1 |
    |     |            |  |     | =  |    |
    |     | 0       S  |  | wk2 |    | x2 |  <<-- x2 already computed.
    |     |            |  |     |    |    |       and we need x1
    |
    |    with x2 = S^{-1} wk2 [assumed to have been computed ] 
    |--------------------------------------------------------------------*/
  /*--------------------  local variables  */
  int j, len=levmat->n, lenB=levmat->nB, *qperm=levmat->perm;
  FLOAT *work = levmat->wk; 
  /*-------------------- copy x onto wk */  
  matvec(levmat->F, &x[lenB], work);   /*  work = F * x_2   */
  arms_Lsol(levmat->L, work, work);         /*  work = L \ work    */
  for (j=0; j<lenB; j++)               /*  wk1 = wk1 - work  */
    work[j] = x[j] - work[j];
  arms_Usol(levmat->U, work, work);         /*  wk1 = U \ wk1 */ 
  memcpy(&work[lenB],&x[lenB],(len-lenB)*sizeof(FLOAT));
  /*---------------------------------------
    |   apply reverse permutation
    |--------------------------------------*/
  for (j=0; j<len; j++)
    wk[j] = work[qperm[j]];
  return 0;
}
/*----end-of-ascend----------------------------------------------------
|----------------------------------------------------------------------
|--------------------------------------------------------------------*/

void matvecz(csptr mata, FLOAT *x, FLOAT *y, FLOAT *z) 
{
  /*---------------------------------------------------------------------
    | This function does the matrix vector  z = y - A x.
    |----------------------------------------------------------------------
    | on entry:
    | mata  = the matrix (in SparRow form)
    | x, y   = two input vector
    |
    | on return
    | z    = the result:  y - A * x
    | z-location must be different from that of x 
    | i.e., y and x are used but not modified.
    |--------------------------------------------------------------------*/
  /*   local variables    */
  int i, k, *ki;
  FLOAT *kr, t;
  for (i=0; i<mata->n; i++) {
    kr = mata->pa[i];
    ki = mata->pj[i];
    t = y[i] ;
    for (k=0; k<mata->nnzrow[i]; k++)
      t -= kr[k] * x[ki[k]];
    z[i] = t; 
  }
  return;
}
/*---------------end of matvecz----------------------------------------
 *--------------------------------------------------------------------*/

p4ptr Lvsol2(FLOAT *x, int nlev, p4ptr levmat, ilutptr ilusch, int
	     flag)  
{
  /* Macro L-solve -- corresponds to left (L) part of arms
     |  preconditioning operation -- 
     |  on entry : 
     |   x =  right- hand side to be operated on by the preconditioner
     |  on return : x is overwritten
     |   x =  output result of operation 
     |  
     |  Note : in-place operation -- b and x can occupy the same space..
     | --------------------------------------------------------------------*/ 
  /*-------------------- local variables  */
  int nloc=levmat->n, first, lenB;
  p4ptr last=levmat; 
  /*-------------------- take care of  special cases :  nlev==0 --> lusol  */
  if (nlev == 0) {
    SchLsol(ilusch,x);
    return (last);
  }
  first = 0;
  /*-------------------- descend                                      */
  while (levmat) { 
    nloc=levmat->n;
    lenB =levmat->nB;
    /*-------------------- left scaling                                  */
    if (levmat->D1 !=  NULL) 
      dscale(nloc,levmat->D1, &x[first],  &x[first]); 
    /*--------------------  RESTRICTION/ DESCENT OPERATION  */
    if (lenB) 
      descend (levmat, &x[first], &x[first]);
    first += lenB; 
    last = levmat;
    levmat = levmat->next;
    /*---------------------------------------------------------------------
      | next level 
      +--------------------------------------------------------------------*/
  }
  if (flag) {
    SchLsol(ilusch,&x[first]);
  }
  return last; 
}

int Uvsol2(FLOAT *x, int nlev, int n, p4ptr levmat,
	   ilutptr ilusch) 
{
  /* Macro U-solve -- corresponds to left (L) part of arms
     |  preconditioning operation -- 
     |  on entry : 
     |  b  =  right- hand side to be operated on by the preconditioner
     |  on return  = x has been overwritten =
     |  x  =  output result of operation 
     |  
     |  Note : in-place operation -- b and x  can occupy the same space..
     | --------------------------------------------------------------------*/ 

  /*-------------------- local variables  */
  int nloc, lenB, first; 
  /*-------------------- work array                                        */
  /*-------------------- take care of  special cases :  nlev==0 --> lusol  */
  /*-------------------- case of zero levels                             */
  if (nlev == 0) { 
    SchUsol(ilusch, x);  
    return(0);
  }
  /*-------------------- general case                               */
  nloc=levmat->n; 
  lenB=levmat->nB; 
  first = n - nloc; 
  /*-------------------- last level                                 */
  first += lenB; 
  SchUsol(ilusch, &x[first]); 
  /*-------------------- other levels                               */
  while (levmat) {
    nloc = levmat->n; 
    first -= levmat->nB;
    if (levmat->n) 
      ascend(levmat, &x[first],&x[first]);
    /*-------------------- right scaling */
    if (levmat->D2 !=  NULL) 
      dscale(nloc, levmat->D2, &x[first], &x[first]) ;
    levmat = levmat->prev; 
  }
  return 0;
  /*--------------------  PROLONGATION/ ASCENT OPERATION */
}

int armsol2(FLOAT *x,  arms Prec) 
{ 
  /* |  combined preconditioning operation -- combines the
     |  left and right actions. 
     | 
     |  on entry : 
     |   x =  right- hand side to be operated on by the preconditioner
     |  on return : x is overwritten - 
     |   x =  output result of operation 
     |  
     |  Note : in-place operation -- b and x can occupy the same space..
     | --------------------------------------------------------------------*/ 
  /*-------------------- local variables  */
  p4ptr levmat = Prec->levmat;
  ilutptr ilusch = Prec->ilus;
  int nlev = Prec->nlev;
  int n = levmat->n; 
  p4ptr last;
  if (nlev == 0) {
    n = ilusch->n;
    SchLsol(ilusch, x); 
    SchUsol(ilusch, x); 
    return 0;
  }
  last = Lvsol2(x, nlev, levmat, ilusch, 1) ;
  Uvsol2(x, nlev, n, last, ilusch) ; 
  return 0; 
}

void SchLsol(ilutptr ilusch, FLOAT *y) 
{
  /*---------------------------------------------------------------------
    |  Forward solve for Schur complement part = 
    |----------------------------------------------------------------------
    | on entry:
    | ilusch  = the LU matrix as provided from the ILU functions.
    | y       = the right-hand-side vector
    |
    | on return
    | y       = solution of LU x = y. [overwritten] 
    |---------------------------------------------------------------------*/
  /*-------------------- local variables                        */
  int n = ilusch->n, j, *perm = ilusch->rperm;
  FLOAT *work = ilusch->wk; 
  /*-------------------- begin: right scaling                          */
  if (ilusch->D1 != NULL) 
    dscale(n, ilusch->D1, y, y); 
  /*-------------------- ONE SIDED ROW PERMS */
  if (perm != NULL) { 
    for (j=0; j<n; j++)
      work[perm[j]] = y[j]; 
    /*--------------------  L solve proper */
    arms_Lsol(ilusch->L, work, y); 
  } else 
    arms_Lsol(ilusch->L, y, y); 
  /*---------------end of SchLsol---------------------------------------
    ----------------------------------------------------------------------*/
}

void SchUsol(ilutptr ilusch, FLOAT *y) 
{
  /*---------------------------------------------------------------------
    | U-solve for Schur complement  - 
    |----------------------------------------------------------------------
    | on entry:
    | ilusch  = the LU matrix as provided from the ILU functions.
    | y       = the right-hand-side vector
    |
    | on return 
    | y       = solution of U x = y. [overwritten on y] 
    |----------------------------------------------------------------------*/
  int n = ilusch->n, j,  *perm = ilusch->perm, *cperm;
  FLOAT *work = ilusch->wk; 
  /* -------------------- begin by U-solving */
  /*-------------------- CASE: column pivoting  used (as in ILUTP) */
  if (ilusch->perm2 != NULL) {
    arms_Usol(ilusch->U, y, y);
    cperm = ilusch->perm2; 
    for (j=0; j<n; j++)
      work[cperm[j]] = y[j];
  }
  else
    /*-------------------- CASE: no column pivoting  used                   */
    arms_Usol(ilusch->U, y, work);
  /*-------------------- generic permutation                              */
  if (perm != NULL) {
    for (j=0; j<n; j++)
      y[j] = work[perm[j]];
  } else
    memcpy(y, work,n*sizeof(FLOAT));
    
  /*-------------------- case when diagonal scaling is done on columns    */ 
  if (ilusch->D2 !=NULL) 
    dscale(n, ilusch->D2, y, y);
}
/*---------------end of SchUsol---------------------------------------
----------------------------------------------------------------------*/

int lusolC( FLOAT *y, FLOAT *x, iluptr lu )
{
  /*----------------------------------------------------------------------
   *    performs a forward followed by a backward solve
   *    for LU matrix as produced by iluk
   *    y  = right-hand-side 
   *    x  = solution on return 
   *    lu = LU matrix as produced by iluk. 
   *--------------------------------------------------------------------*/
  int n = lu->n, i, j, nnzrow, *ja;
  FLOAT *D;
  csptr L, U;

  L = lu->L;
  U = lu->U;
  D = lu->D;

  /* Block L solve */
  for( i = 0; i < n; i++ ) {
    x[i] = y[i];
    nnzrow = L->nnzrow[i];
    ja = L->pj[i];
    for( j = 0; j < nnzrow; j++ ) {
      x[i] -= x[ja[j]] * L->pa[i][j];
    }
  }
  /* Block -- U solve */
  for( i = n-1; i >= 0; i-- ) {
    nnzrow = U->nnzrow[i];
    ja = U->pj[i];
    for( j = 0; j < nnzrow; j++ ) {
      x[i] -= x[ja[j]] * U->pa[i][j];
    }
    x[i] *= D[i];
  }
  return (0); 
}

int rpermC(csptr mat, int *perm)
{
  /*----------------------------------------------------------------------
    |
    | This subroutine permutes the rows of a matrix in SparRow format. 
    | rperm  computes B = P A  where P is a permutation matrix.  
    | The permutation P is defined through the array perm: for each j, 
    | perm[j] represents the destination row number of row number j. 
    |
    |-----------------------------------------------------------------------
    | on entry:
    |----------
    | (amat) = a matrix stored in SparRow format.
    |
    |
    | on return:
    | ----------
    | (amat) = P A stored in SparRow format.
    |
    | integer value returned:
    |             0   --> successful return.
    |             1   --> memory allocation error.
    |---------------------------------------------------------------------*/
  int **addj, *nnz, i, size=mat->n;
  FLOAT **addm;
  addj = (int **)Malloc( size*sizeof(int *), "rpermC" );
  addm = (FLOAT **) Malloc( size*sizeof(FLOAT *), "rpermC" );
  nnz = (int *) Malloc( size*sizeof(int), "rpermC" );
  for (i=0; i<size; i++) {
    addj[perm[i]] = mat->pj[i];
    addm[perm[i]] = mat->pa[i];
    nnz[perm[i]] = mat->nnzrow[i];
  }
  for (i=0; i<size; i++) {
    mat->pj[i] = addj[i];
    mat->pa[i] = addm[i];
    mat->nnzrow[i] = nnz[i];
  }
  free(addj);
  free(addm);
  free(nnz);
  return 0;
}

int cpermC(csptr mat, int *perm) 
{
  /*----------------------------------------------------------------------
    |
    | This subroutine permutes the columns of a matrix in SparRow format.
    | cperm computes B = A P, where P is a permutation matrix.
    | that maps column j into column perm(j), i.e., on return 
    | The permutation P is defined through the array perm: for each j, 
    | perm[j] represents the destination column number of column number j. 
    |
    |-----------------------------------------------------------------------
    | on entry:
    |----------
    | (mat) = a matrix stored in SparRow format.
    |
    |
    | on return:
    | ----------
    | (mat) = A P stored in SparRow format.
    |
    | integer value returned:
    |             0   --> successful return.
    |             1   --> memory allocation error.
    |---------------------------------------------------------------------*/
  int i, j, *newj, size=mat->n, *aja;
  newj = (int *) Malloc( size*sizeof(int), "cpermC" );
  for (i=0; i<size; i++) {
    aja = mat->pj[i];
    for (j=0; j<mat->nnzrow[i]; j++)
      newj[j] = perm[aja[j]];
  
    for (j=0; j<mat->nnzrow[i]; j++)
      aja[j] = newj[j];
    mat->pj[i] = aja;
  }
  free(newj);
  return 0;
}

int dpermC(csptr mat, int *perm) 
{
  /*----------------------------------------------------------------------
    |
    | This subroutine permutes the rows and columns of a matrix in 
    | SparRow format.  dperm computes B = P^T A P, where P is a permutation 
    | matrix.
    |
    |-----------------------------------------------------------------------
    | on entry:
    |----------
    | (amat) = a matrix stored in SparRow format.
    |
    |
    | on return:
    | ----------
    | (amat) = P^T A P stored in SparRow format.
    |
    | integer value returned:
    |             0   --> successful return.
    |             1   --> memory allocation error.
    |---------------------------------------------------------------------*/
  if (rpermC(mat, perm)) return 1;
  if (cpermC(mat, perm)) return 1;
  return 0;
}

int CSparTran( csptr amat, csptr bmat, CompressType *compress )
{
  /*----------------------------------------------------------------------
    | Finds the compressed transpose of a matrix stored in SparRow format.
    | Patterns only.
    |-----------------------------------------------------------------------
    | on entry:
    |----------
    | (amat)     = a matrix stored in SparRow format.
    | (compress) = quotient graph of matrix amat
    |
    | on return:
    | ----------
    | (bmat)     = the compressed transpose of (mata) stored in SparRow
    |              format.
    |
    | integer value returned:
    |             0   --> successful return.
    |             1   --> memory allocation error.
    |---------------------------------------------------------------------*/
  int i, j, *ind, nnzrow, pos, size=amat->n, *aja;
  ind = bmat->nnzrow;

  for (i=0; i<size; i++)
    ind[i] = 0;
  /*-------------------- compute lengths  */
  for (i=0; i<size; i++) {
    if( compress[i].grp != -1 ) continue;
    aja = amat->pj[i];
    nnzrow = amat->nnzrow[i];
    for (j=0; j < nnzrow; j++) {
      pos = aja[j];
      if( compress[pos].grp == -1 ) {
	ind[pos]++;
      }
    }
  }

  /*--------------------  allocate space  */
  for (i=0; i<size; i++) {
    if( ind[i] == 0 ) {
      bmat->pj[i] = NULL;
      continue;
    }
    bmat->pj[i] = (int *)Malloc( ind[i]*sizeof(int), "CSparTran" );
    ind[i] = 0; /* indicate next available position of each row */
  }
  /*--------------------  now do the actual copying  */
  for (i=0; i<size; i++) {
    if( compress[i].grp != -1 ) continue;
    aja = amat->pj[i];
    nnzrow = amat->nnzrow[i];
    for (j = 0; j < nnzrow; j++) {
      pos = aja[j];
      if( compress[pos].grp == -1 ) {
	bmat->pj[pos][ind[pos]] = i;
	ind[pos]++;
      }
    }
  }
  return 0;
}

int condestArms(arms armspre, FLOAT *y, FILE *fp )
{
  /*-------------------- simple estimate of cond. number of precon */
  int n = armspre->n, i;
  double norm = 0.0;
  
  for( i = 0; i < n; i++ ) 
    y[i] = 1.0;
  armsol2(y, armspre)  ;
  for( i = 0; i < n; i++ ) {
    norm = max( norm, ABS_VALUE(y[i]) );
  }
  fprintf( fp, "ARMS inf-norm lower bound = : %16.2f\n", norm );
  if( norm > 1e30 ) {
    return -1;
  }
  return 0;
}

void Lsolp(int start, csptr mata, FLOAT *b, FLOAT *y)
{
  /*---------------------------------------------------------------------
    |
    | This routine does the forward solve L y = b. where L is upper or
    | bottom part of local low triangular matrix
    |
    | Can be done in place.
    |
    | Zhongze Li, Aug. 17th, 2001
    |
    |----------------------------------------------------------------------
    | on entry:
    | start = the index of the first component
    | n     = one ore than the index owned by the processor 
    | b     = a vector
    | mata  = the matrix (in SparRow form)
    |
    | on return
    | y     = the product L^{-1} * b
    |
    |--------------------------------------------------------------------*/
  /*   local variables    */
  int i, k;
  FLOAT *kr;
  int *ki, n;

  n = mata->n;
  for (i=0; i<n; i++) {
    y[i] = b[i];
    if ( mata->nnzrow[i] > 0 ) {
      kr = mata->pa[i];
      ki = mata->pj[i];
      for (k=0; k<mata->nnzrow[i]; k++)
	if(ki[k] < start) {
	  y[i] -= kr[k]*y[ki[k]];
	}
    }
  }
}

void Usolp(int start, csptr mata, FLOAT *y, FLOAT *x)
{
  /*---------------------------------------------------------------------
    |
    | This routine does the backward solve U x = y, where U is upper or
    | bottom part of local upper triangular matrix
    |
    | Can be done in place.
    |
    | Zhongze Li, Aug. 17th, 2001
    |
    |----------------------------------------------------------------------
    | on entry:
    | start = the index of the first component
    | n     = one ore than the index owned by the processor 
    | y     = a vector
    | mata  = the matrix (in SparRow form)
    |
    | on return
    | x     = the product U^{-1} * y
    |
    |---------------------------------------------------------------------*/
  /*   local variables    */
  int i, k, *ki;
  FLOAT *kr, t;

  for (i=start-1; i>= 0; i--) {
    kr = mata->pa[i];
    ki = mata->pj[i];
    t  = y[i];
    for (k=1; k<mata->nnzrow[i]; k++)
      t -= kr[k] * x[ki[k]];
    x[i] = t*kr[0];
  }
}
