#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#include "protos.h"

void parms_errexit( char *f_str, ... )
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s\n", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  exit( -1 );
}

void *Malloc( int nbytes, char *msg )
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL)
    parms_errexit( "Not enough mem for %s. Requested size: %d bytes", msg, nbytes );

  return ptr;
}

int setupCS(csptr amat, int len, int job)
{
/*----------------------------------------------------------------------
| Initialize SparRow structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SparRow struct.
|     len   =  size of matrix
|     job   =  0: pattern only
|              1: data and pattern
|
| On return:
|===========
|
|  amat->n
|      ->*nnzrow
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   amat->n = len;
   amat->nnzrow = (int *)Malloc( len*sizeof(int), "setupCS" );
   amat->pj = (int **) Malloc( len*sizeof(int *), "setupCS" );
   if( job == 1 ) 
       amat->pa = (FLOAT **) Malloc( len*sizeof(FLOAT *), "setupCS" );
   else
       amat->pa = NULL;
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupCS
|--------------------------------------------------------------------*/

int cleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SparRow structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SparRow struct.
|--------------------------------------------------------------------*/
   /*   */
  int i;
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;

  for (i=0; i<amat->n; i++) {
    if (amat->nnzrow[i] > 0) {
      if( amat->pa ) free(amat->pa[i]);
      free(amat->pj[i]);
    }
  }    
  if (amat->pa) free(amat->pa);
  free(amat->pj);
  free(amat->nnzrow);
  free(amat);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanCS
|--------------------------------------------------------------------*/

int nnzCS( csptr amat )
{
    int nnz = 0, i, n = amat->n;
    for( i = 0; i < n; i++ ) {
        nnz += amat->nnzrow[i];
    }
    return nnz;
}

int cscpy(csptr amat, csptr bmat){
/*----------------------------------------------------------------------
| Convert CSR matrix to SparRow struct
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )   = Matrix stored in SparRow format
|
|
| On return:
|===========
|
| ( bmat )  =  Matrix stored as SparRow struct containing a copy
|              of amat 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int j, len, size=amat->n;
  FLOAT *bma;
  int *bja;
/*------------------------------------------------------------*/
  for (j=0; j<size; j++) {
    len = bmat->nnzrow[j] = amat->nnzrow[j];
    if (len > 0) {
      bja = (int *) Malloc(len*sizeof(int), "cscpy:1" );
      bma = (FLOAT *) Malloc(len*sizeof(FLOAT), "cscpy:2" );
      memcpy(bja,amat->pj[j],len*sizeof(int));
      memcpy(bma,amat->pa[j],len*sizeof(FLOAT));
      bmat->pj[j] = bja;
      bmat->pa[j] = bma;
    }	
  }
  return 0;
}
/*-----------------------------------------------------------------------*/

int setupILU( iluptr lu, int n )
{
/*----------------------------------------------------------------------
| Initialize ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|       n   =  size of matrix
|
| On return:
|===========
|
|    lu->n
|      ->L     L matrix, SparRow format
|      ->D     Diagonals
|      ->U     U matrix, SparRow format
|      ->work  working buffer of length n
|      ->bf    buffer
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    lu->n  = n;
    lu->D = (FLOAT *)Malloc( sizeof(FLOAT) * n, "setupILU" );
    lu->L = (csptr)Malloc( sizeof(SparMat), "setupILU" );
    setupCS( lu->L, n, 1 );
    lu->U = (csptr)Malloc( sizeof(SparMat), "setupILU" );
    setupCS( lu->U, n, 1 );
    lu->work = (int *)Malloc( sizeof(int) * n, "setupILU" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of setupILU
|--------------------------------------------------------------------*/


int cleanILU( iluptr lu )
{
/*----------------------------------------------------------------------
| Free up memory allocated for ILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|--------------------------------------------------------------------*/
  if( NULL == lu ) return 0;
  if( lu->D ) {
    free( lu->D );
  }
  cleanCS( lu->L );
  cleanCS( lu->U );
  if( lu->work ) free( lu->work );
  free( lu );
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILU
|--------------------------------------------------------------------*/

int mallocRow( iluptr lu, int nrow )
{
/*----------------------------------------------------------------------
| Prepare space of a row according to the result of level structure
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a ILUSpar struct.
|     nrow  =  the current row to deal with
|
| On return:
|===========
|
|    lu->L->pa[nrow][...]
|      ->U->pa[nrow][...]
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int nnzrow = lu->L->nnzrow[nrow];
    lu->L->pa[nrow] = (FLOAT *)Malloc( sizeof(FLOAT)*nnzrow, "mallocRow" );
    nnzrow = lu->U->nnzrow[nrow];
    lu->U->pa[nrow] = (FLOAT *)Malloc( sizeof(FLOAT)*nnzrow, "mallocRow" );
    return 0;
}
/*---------------------------------------------------------------------
|     end of mallocRow
|--------------------------------------------------------------------*/

int setupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E) 
{
/*----------------------------------------------------------------------
| initialize PerMat4 struct given the F, E, blocks.  
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a PerMat4 struct.
|     Bn    =  size of B block
|     Cn    =  size of C block
|     F, E  = the two blocks to be assigned to srtruct - without the
|
| On return:
|===========
|
|  amat->L                for each block: amat->M->n
|      ->U                                       ->nnzrow
|      ->E                                       ->pj
|      ->F                                       ->pa
|      ->perm
|      ->rperm       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int n;
   /* size n */
   n = amat->n = Bn + Cn;
   amat->nB = Bn; 
/* amat->perm = (int *) Malloc(n*sizeof(int), "setupP4:1" ); */
/*   assign space for wk -- note that this is only done at 1st level
     at other levels, copy pointer of wk from previous level */
   if (amat->prev == NULL)  /* wk has 2 * n entries now */
     amat->wk   = (FLOAT *) Malloc(2*n*sizeof(FLOAT), "setupP4:2" );
   else 
     amat->wk = (amat->prev)->wk; 

/*-------------------- L and U */ 
   amat->L = (csptr) Malloc(sizeof(SparMat), "setupP4:3" );
   if (setupCS(amat->L, Bn,1)) return 1;
   /*    fprintf(stdout,"  -- BN %d   Cn   %d \n", Bn,Cn);  */
   amat->U = (csptr) Malloc(sizeof(SparMat), "setupP4:4" );
   if (setupCS(amat->U, Bn,1)) return 1;

   amat->F = F; 
   amat->E = E; 
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupP4 
|--------------------------------------------------------------------*/

int cleanP4(p4ptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for Per4Mat structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
|--------------------------------------------------------------------*/

/*  -------------------------- */
  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
   

  if (amat->perm) {
    free(amat->perm); 
    amat->perm = NULL;
  }
  
  if (!amat->symperm) { 
    if (amat->rperm) free(amat->rperm); 
    amat->rperm = NULL;
  } 
  
  if (amat->F) {
    cleanCS(amat->F); 
    amat->F = NULL;
  }
  if (amat->E) {
    cleanCS(amat->E); 
    amat->E = NULL;
  }
  if (amat->L) {
    cleanCS(amat->L);
    amat->L = NULL;
   }
  if (amat->U) {
    cleanCS(amat->U);
    amat->U = NULL;
  }
  
  if (amat->prev == NULL) 
    if (amat->wk) free(amat->wk);  
  
  if (amat->D1) free(amat->D1);
  if (amat->D2) free(amat->D2);
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanP4
|--------------------------------------------------------------------*/


int setupILUT(ilutptr amat, int len)
{
/*----------------------------------------------------------------------
| Allocate pointers for ILUTfac structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a ILUTfac struct.
|     len   =  size of L U  blocks
|
| On return:
|===========
|
|  amat->L                for each block: amat->M->n
|      ->U                                       ->nnzrow
|                                                ->pj
|                                                ->pa
|      ->rperm       (if meth[0] > 0)
|      ->perm2       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Permutation arrays are initialized to the identity.
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  amat->n = len;
 amat->wk = (FLOAT *) Malloc(2*len*sizeof(FLOAT), "setupILUT:5" );
  amat->L = (csptr) Malloc(sizeof(SparMat), "setupILUT:6" );
  if (setupCS(amat->L, len,1)) return 1;
  amat->U = (csptr) Malloc(sizeof(SparMat), "setupILUT:7" );
  if (setupCS(amat->U, len,1)) return 1;
  return 0;
}    
/*---------------------------------------------------------------------
|     end of setupILUT
|--------------------------------------------------------------------*/
int cleanILUT(ilutptr amat, int indic)
{
/*----------------------------------------------------------------------
| Free up memory allocated for IluSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a IluSpar struct.
|  indic    = indicator for number of levels.  indic=0 -> zero level.
|--------------------------------------------------------------------*/
  /*----------------*/
   
  if (amat->wk) {
    free(amat->wk); 
    amat->wk = NULL;
  }
  cleanCS(amat->L);
  cleanCS(amat->U);

  if (indic) cleanCS(amat->C);  
/*-------------------- nonsymmetric permutation */
  if (amat->rperm) {
    free(amat->rperm);
    amat->rperm = NULL;
  }
  if (amat->perm) {
    free(amat->perm); 
    amat->perm = NULL;
  }
  
/*-------------------- ilutp permutation */
  if (amat->perm2) free(amat->perm2);
/*-------------------- diagonal scalings */
   if (amat->D1) free(amat->D1);
   if (amat->D2) free(amat->D2);
   return 0;
}
/*---------------------------------------------------------------------
|     end of cleanILUT
|--------------------------------------------------------------------*/


void setup_arms (arms Levmat) 
{
  Levmat->ilus = (ilutptr) Malloc(sizeof(IluSpar), "setup_arms:ilus" );
  Levmat->levmat = (p4ptr) Malloc(sizeof(Per4Mat), "setup_arms:levmat" );
}
int cleanARMS(arms ArmsPre)
{
  p4ptr amat = ArmsPre->levmat;
  ilutptr cmat = ArmsPre->ilus;
/*----------------------------------------------------------------------
| Free up memory allocated for entire ARMS preconditioner.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
| ( cmat )  =  Pointer to a IluSpar struct.
|--------------------------------------------------------------------*/
/* case when nlev == 0 */  
  int indic=(amat->nB != 0) ;
    /*  && amat->next !=NULL) ; */
  
  p4ptr levc, levn;
/*
  int cleanCS(csptr);
  int cleanP4(p4ptr);
  int cleanILUT(ilutptr, int);
*/
  levc = amat; 

  if (indic) { 
    while (levc) {
      if (cleanP4(levc)) return(1) ; 
      levn = levc->next;
      free(levc);
      levc = levn;
    }		
  }	
   else 	
     if (amat) {
       free(amat) ; 
       amat = NULL;
     }
  
  cleanILUT(cmat,indic); 
  
  
  if (cmat) {
    free(cmat);	
    cmat = NULL;
  }
  
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanARMS 
|--------------------------------------------------------------------*/


int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
	     csptr E, csptr C)
{
/*---------------------------------------------------------------------
| Convert permuted csrmat struct to PerMat4 struct 
|                - matrix already permuted
|----------------------------------------------------------------------
| on entry:
|========== 
| ( amat )  =  Matrix stored in SparRow format.
|              Internal pointers (and associated memory) destroyed before
|              return.
|
| On return:
|===========
|
| B, E, F, C = 4 blocks in 
| 
|          | B   F |      
|   Amat = |       | 
|          | E   C | 
| 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int j, j1, numr, numl, ind, newj, rowz, *rowj, *new1j, *new2j;
   FLOAT *rowm, *new1m, *new2m;
/*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/
   if (setupCS(B,bsize,1)) goto label111; 
   if (setupCS(F,bsize,1)) goto label111;
   if (setupCS(E,csize,1)) goto label111;
   if (setupCS(C,csize,1)) goto label111;
   new1j = (int *) Malloc(bsize*sizeof(int), "csSplit4:1" );
   new2j = (int *) Malloc(csize*sizeof(int), "csSplit4:2" );
   new1m = (FLOAT *) Malloc(bsize*sizeof(FLOAT), "csSplit4:3" );
   new2m = (FLOAT *) Malloc(csize*sizeof(FLOAT), "csSplit4:4" );
/*    B and F blocks */ 
   for (j=0; j<bsize; j++) {
      numl = numr = 0;
      rowz = amat->nnzrow[j];
      rowj = amat->pj[j];
      rowm = amat->pa[j];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      B->nnzrow[j] = numl;
      F->nnzrow[j] = numr;
      if (numl>0) {
	 B->pj[j] = (int *) Malloc(numl*sizeof(int), "csSplit4:5" );
	 B->pa[j] = (FLOAT *) Malloc(numl*sizeof(FLOAT), "csSplit4:6" );
      }
      if (numr>0) {
	 F->pj[j] = (int *) Malloc(numr*sizeof(int), "csSplit4:7" );
	 F->pa[j] = (FLOAT *) Malloc(numr*sizeof(FLOAT), "csSplit4:8" );
      }
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	 newj = rowj[j1];
	 if (newj<bsize) {
	    new1j[numl] = newj;
	    new1m[numl] = rowm[j1];
	    numl++;
	 }
	 else {
	    new2j[numr] = newj - bsize;
	    new2m[numr] = rowm[j1];
	    numr++;
	 }
      }
      memcpy(B->pj[j], new1j, numl*sizeof(int));
      memcpy(B->pa[j], new1m, numl*sizeof(FLOAT));
      memcpy(F->pj[j], new2j, numr*sizeof(int));
      memcpy(F->pa[j], new2m, numr*sizeof(FLOAT));
   }
/*    E and C blocks */
   for (j=0; j<csize; j++) {
      numl = numr = 0;
      ind = bsize + j;
      rowz = amat->nnzrow[ind];
      rowj = amat->pj[ind];
      rowm = amat->pa[ind];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      E->nnzrow[j] = numl;
      C->nnzrow[j] = numr;
      if (numl>0) {
	E->pj[j] = (int *) Malloc(numl*sizeof(int), "csSplit4:9" );
	E->pa[j] = (FLOAT *) Malloc(numl*sizeof(FLOAT), "csSplit4:10" );
      }	
      if (numr>0) {
	C->pj[j] = (int *) Malloc(numr*sizeof(int), "csSplit4:11" );
	C->pa[j] = (FLOAT *) Malloc(numr*sizeof(FLOAT), "csSplit4:12" );
      }		
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	newj = rowj[j1];
	if (newj<bsize) {
	  new1j[numl] = newj;
	  new1m[numl] = rowm[j1];
	  numl++;
	}
	else {
	  new2j[numr] = newj - bsize;
	  new2m[numr] = rowm[j1];
	  numr++;
	}
      }
      memcpy(E->pj[j], new1j, numl*sizeof(int));
      memcpy(E->pa[j], new1m, numl*sizeof(FLOAT));
      memcpy(C->pj[j], new2j, numr*sizeof(int));
      memcpy(C->pa[j], new2m, numr*sizeof(FLOAT));
   }

   if (new1j) free(new1j);
   if (new2j) free(new2j);
   if (new1m) free(new1m);
   if (new2m) free(new2m);
   return 0;
label111:
   return 1;
}
/*---------------------------------------------------------------------
|     end of csSplit4
|--------------------------------------------------------------------*/

int CSRcs( int n, FLOAT *a, int *ja, int *ia, csptr mat, int rsa )
{
/*----------------------------------------------------------------------
| Convert CSR matrix to SparRow struct
|----------------------------------------------------------------------
| on entry:
|==========
| a, ja, ia  = Matrix stored in CSR format (with FORTRAN indexing).
| rsa        = source file is symmetric HB matrix 
|
| On return:
|===========
|
| ( mat )  =  Matrix stored as SparRow struct.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int i, j, j1, len, col, nnz;
  FLOAT *bra;
  int *bja;
  /*    setup data structure for mat (csptr) struct */
  setupCS( mat, n, 1 );

  if( rsa ) { /* RSA HB matrix */
    for( j = 0; j < n; j++ ) {
      len = ia[j+1] - ia[j];
      mat->nnzrow[j] = len;
    }
    for( j = 0; j < n; j++ ) {
      for( j1 = ia[j]-1; j1 < ia[j+1]-1; j1++ ) {
        col = ja[j1] - 1;
        if( col != j ) mat->nnzrow[col]++;
      }
    }
    for( j = 0; j < n; j++ ) {
      nnz = mat->nnzrow[j];
      mat->pj[j] = (int *)Malloc( nnz * sizeof(int), "CSRcs" );
      mat->pa[j] = (FLOAT *)Malloc( nnz * sizeof(FLOAT), "CSRcs" );
      mat->nnzrow[j] = 0;
    }
    for( j = 0; j < n; j++ ) {
      for( j1 = ia[j]-1; j1 < ia[j+1]-1; j1++ ) {
        col = ja[j1] - 1;
        mat->pj[j][mat->nnzrow[j]] = col;
        mat->pa[j][mat->nnzrow[j]] = a[j1];
        mat->nnzrow[j]++;
        if( col != j ) {
          mat->pj[col][mat->nnzrow[col]] = j;
          mat->pa[col][mat->nnzrow[col]] = a[j1];
          mat->nnzrow[col]++;
        }
      }
    }
    return 0;
  }
  for (j=0; j<n; j++) {
    len = ia[j+1] - ia[j];
    mat->nnzrow[j] = len;
    if (len > 0) {
      bja = (int *) Malloc( len*sizeof(int), "CSRcs" );
      bra = (FLOAT *) Malloc( len*sizeof(FLOAT), "CSRcs" );
      i = 0;
      for (j1=ia[j]-1; j1<ia[j+1]-1; j1++) {
        bja[i] = ja[j1] - 1;
        bra[i] = a[j1] ;
        i++;
      }
      mat->pj[j] = bja;
      mat->pa[j] = bra;
    }
  }    
  return 0;
}
/*---------------------------------------------------------------------
|     end of CSRcs
|--------------------------------------------------------------------*/

/*
int nnz_ilu( iluptr lu )
{
  int nnz = 0, i;
  for( i = 0; i < lu->n; i++ ) {
    nnz += lu->L->nnzrow[i];
    nnz += lu->U->nnzrow[i];
    nnz++;
  }
  return nnz;
}
*/

int lev4_nnz(p4ptr levmat, int *lev) 
{
  /* counts all nonzero elements in levmat struct  -- 
     recursive */
  int nnzT, nnzL, nnzU, nnzF, nnzE, nnzDown=0;
  p4ptr nextmat; 
  nnzL = cs_nnz(levmat->L); 
  nnzU = cs_nnz(levmat->U); 
  nnzF = cs_nnz(levmat->F); 
  nnzE = cs_nnz(levmat->E); 
  nnzT = nnzL+nnzU+nnzF+nnzE;
  /* print */
#if 0
  if (*lev == 0) 
    fprintf(ft,
	    "\nnnz/lev used:      L        U        F        E    subtot\n");  
  fprintf(ft,"    Level %2d %8d %8d %8d %8d %8d\n",
	  *lev, nnzL, nnzU, nnzF, nnzE, nnzT);
#endif 
  (*lev)++;
  nextmat = levmat->next; 
  if (nextmat != NULL) 
   nnzDown = lev4_nnz(nextmat, lev);
   return (nnzT+nnzDown); 
}
  
int cs_nnz (csptr A) 
{
  /* counts the number of nonzero elements in CSR matrix A */
  int i, n, nnz=0; 
  n = A->n; 
  for (i=0; i<n; i++) nnz +=A->nnzrow[i];
  return nnz;
  }

int nnz_arms (arms PreSt)
{ 
/*-------------------------------------------------------
| computes and prints out total number of nonzero elements
| used in ARMS factorization 
+--------------------------------------------------------*/
  p4ptr levmat   = PreSt->levmat; 
  ilutptr ilschu = PreSt->ilus; 
  int nlev       = PreSt->nlev;
  int ilev=0,nnz_lev,nnz_sch,nnz_tot; 
  nnz_lev = 0; 
  if (nlev) nnz_lev+= lev4_nnz(levmat, &ilev);
  nnz_sch = cs_nnz(ilschu->L)+cs_nnz(ilschu->U);
  if (nlev) nnz_sch += cs_nnz(ilschu->C);
  nnz_tot = nnz_lev+nnz_sch;
#if 0
  fprintf(ft,"\n"); 
  fprintf(ft,"Total nonzeros for interm. blocks.... =  %10d\n",nnz_lev);
  fprintf(ft,"Total nonzeros for last level ....... =  %10d\n",nnz_sch);
  fprintf(ft,"Grand total.......................... =  %10d\n",nnz_tot);
#endif 
  return nnz_tot;
}

/*----------------------------------------------------------------------
| Output the pattern of L\U, which can be loaded by matlab
----------------------------------------------------------------------*/
/*
int outputLU( iluptr lu, char *filename )
{
  FILE *fmatlab = fopen( filename, "w" );
  int n = lu->n, i, j, nnzrow;
  csptr L = lu->L, U = lu->U;
  
  if( !fmatlab ) return -1;
  fprintf( fmatlab, "%d %d 0\n", n, n );
  for( i = 0; i < n; i++ ) {
    nnzrow = L->nnzrow[i];
    for( j = 0; j < nnzrow; j++ )
      fprintf( fmatlab, "%d %d 1\n", i+1, L->pj[i][j]+1 );
  }
  for( i = 0; i < n; i++ ) {
    nnzrow = U->nnzrow[i];
    for( j = 0; j < nnzrow; j++ )
      fprintf( fmatlab, "%d %d 1\n", i+1, U->pj[i][j]+1 );
  }
  for( i = 0; i < n; i++ )
    fprintf( fmatlab, "%d %d 1\n", i+1, i+1 );
  fclose( fmatlab );
  return 0;
}
*/

