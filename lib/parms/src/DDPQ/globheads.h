#ifndef __VBLOCK_HEADER_H__
#define __VBLOCK_HEADER_H__

#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#include "parms_mat_impl.h"
#include "parms_opt_impl.h"

#define MAX_BLOCK_SIZE   100

/* FORTRAN style vblock format, compatible for many FORTRAN routines */
#define DATA(a,row,i,j)  (a[(j)*(row)+(i)])

/* the dimension of ith Block */
#define B_DIM(bs,i)      (bs[i+1]-bs[i])

#if 0
typedef struct SparRow {
/*--------------------------------------------- 
| C-style CSR format - used internally
| for all matrices in CSR format 
|---------------------------------------------*/
  int n;
  int *nnzrow;  /* length of each row */
  int **ja;     /* pointer-to-pointer to store column indices  */
  double **ma;  /* pointer-to-pointer to store nonzero entries */
} SparMat, *csptr;
#endif 

typedef struct parms_vcsr SparMat;
typedef parms_vcsr csptr;

typedef struct ILUfac {
    int n;
    csptr L;      /* L part elements                            */
    FLOAT *D;    /* diagonal elements                          */
    csptr U;      /* U part elements                            */
    int *work;    /* working buffer */
} ILUSpar, LDUmat, *iluptr;


typedef struct PerMat4 *p4ptr;
typedef struct PerMat4 {
/*------------------------------------------------------------
| struct for storing the block LU factorization 
| contains all the block factors except the 
| data related to the last block. 
| n       = size of current block
| symperm = whether or not permutations are symmetric.
|           used only in cleanP4..
| nB      = size of B-block
| L, U    = ILU factors of B-block
| F, E    = sparse matrices in (1,2) and (2,1) 
|           parts of matrix. 
| perm    = (symmetric) permutation used in factorization
|           comes from the independent set ordering
| rperm   = unsymmetric permutation (rows) not used in this
|           version -- but left here for compatibility..
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|----------------------------------------------------------*/ 
  int n;                  
  int nB; 
  int symperm;
/*   LU factors  */
  csptr L;
  csptr U;
/* E, F blocks   */
  csptr E;
  csptr F;
  int *rperm;       /* row permutation         */ 
  int *perm;        /* col. permutation        */ 
  double *D1 ;      /* diagonal scaling row    */  
  double *D2 ;      /* diagonal scaling columns*/  
  FLOAT *wk;       /* work array              */
/* pointer to next and previous struct         */
  p4ptr prev; 
  p4ptr next;
} Per4Mat; 
/* -------------------------------------------------------------------*/
typedef struct ILUTfac *ilutptr;
typedef struct ILUTfac {
/*------------------------------------------------------------
| struct for storing data related to the last schur complement 
| we need to store the C matrix associated with the last block
| and the ILUT factorization of the related Schur complement.
| 
| n       = size of C block = size of Schur complement
| C       = C block of last level matrix. 
| L, U    = ILU factors of last schur complement. 
|
| meth[4] = parameters for defining variants in factorization 
|           - see function readin for details
| rperm    = row permutation used for very nonsymmetric matrices 
|            [such as bottleneck transversal] -- NOT IN THIS VERSION
| perm2     = unsymmetric permutation (columns) - used primarily
|           for the ILUTP version of ILUT/.. 
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|-----------------------------------------------------------*/
   int n;                  
 /*-------------------- C matrix of last block */
   csptr C;
  /* LU factorization       */
   csptr L;
   csptr U;
 /*--------------------  related to variants and methods */
 /*    int meth[4];   */
  int *rperm;   /* row single-sinded permutation */
  int *perm;    /* column perm .                */
  int *perm2;   /* column permutation coming from pivoting in ILU */ 
  double *D1;
  double *D2;
  FLOAT *wk;
} IluSpar;

typedef struct parms_arms_data {
  int n;			/* dimension of matrix */
  int nlev;			/* number of levels */
  ilutptr ilus;
  p4ptr levmat;
  int ipar[18];
  double pgfpar[2];
  int schur_start;
  int ind;
  int nnz_mat;
  int nnz_prec;
} *parms_arms_data;

typedef struct parms_arms_data armsMat;
typedef parms_arms_data arms;

typedef struct __CompressType
{
  int grp;   /* -1: begin new group, >=0: belong to grp-th row */
  int count; /* block size, valid only if grp = -1 */
} CompressType;

#endif  /* __VBLOCK_HEADER_H__ */

