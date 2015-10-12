/*----------------------------------------------------------------------
 * Parallel Multi-Level Block ILUT PRECONDITIONER
 * arms2    :  parallel arms2
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#define  PERMTOL  0.99   /*  0 --> no permutation 0.01 to 0.1 good  */
#include "protos.h" 

static int arms_free_vcsr(parms_Operator *self)
{
  parms_arms_data data;

  data = (parms_arms_data)(*self)->data;
  cleanARMS(data);
  PARMS_FREE(data);
  return 0;
}

static int arms_view_vcsr(parms_Operator self, parms_Viewer v)
{
/* So far, only available for viewing last level operator */

  parms_arms_data data;
  ilutptr LU;
  int i, j, n, nnz, *pj;
  FLOAT *pa;
  FILE *fp;

  parms_ViewerGetFP(v, &fp);
  data = (parms_arms_data)self->data;
  LU = data->ilus;
  n = LU->n;

/* L part */
  fprintf(fp, "L part of the last level matrix:\n");
  fprintf(fp, "n = %d\n", n);
#if defined(DBL_CMPLX)  
  for (i = 0; i < n; i++) {
    nnz = LU->L->nnzrow[i];
    pj  = LU->L->pj[i];
    pa  = LU->L->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f  %f) ", i, pj[j], creal(pa[j]), cimag(pa[j]));
    }
  }
#else
  for (i = 0; i < n; i++) {
    nnz = LU->L->nnzrow[i];
    pj  = LU->L->pj[i];
    pa  = LU->L->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f) ", i, pj[j], pa[j]);
    }
  }
#endif
/* U part */
  fprintf(fp, "U part of the matrix:\n");
  fprintf(fp, "n = %d\n", n);
#if defined(DBL_CMPLX)
  for (i = 0; i < n; i++) {
    nnz = LU->U->nnzrow[i];
    pj  = LU->U->pj[i];
    pa  = LU->U->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f  %f) ", i, pj[j], creal(pa[j]), cimag(pa[j]));
    }
  }
#else
  for (i = 0; i < n; i++) {
    nnz = LU->U->nnzrow[i];
    pj  = LU->U->pj[i];
    pa  = LU->U->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f) ", i, pj[j], pa[j]);
    }
  }  
#endif
  parms_ViewerStoreFP(v, fp);

  return 0;
}

static void parms_arms_nnz(parms_Operator self, int *nnz_mat, int
			   *nnz_pc)
{
  parms_arms_data data;

  data = (parms_arms_data)self->data;
  *nnz_mat = data->nnz_mat;
  *nnz_pc  = data->nnz_prec;
}

static int parms_arms_lsol_vcsr(parms_Operator self, FLOAT *y, FLOAT
				*x)
{
  parms_arms_data data;
  p4ptr levmat;
  ilutptr ilus;
  int *ipar, schur_start, nlev, n;

  data = (parms_arms_data)self->data;
  levmat = data->levmat;
  ilus   = data->ilus;
  ipar   = data->ipar;
  schur_start = data->schur_start;

  if (ipar[0] == 0) {
    Lsolp(schur_start, ilus->L, y, x); 
    return 0;
  }
  nlev = data->nlev;
  n = data->n;
  PARMS_MEMCPY(x, y, n);
  Lvsol2(x, nlev, levmat, ilus, 0) ;
  return 0;
}

static int parms_arms_sol_vcsr(parms_Operator self, FLOAT *y,
			  FLOAT *x) 
{
  parms_arms_data data;
  int n;

  data = (parms_arms_data)self->data;

  n = data->n;

  PARMS_MEMCPY(x, y, n);

  armsol2(x, data);

  return 0;
}

static int parms_arms_invs_vcsr(parms_Operator self,  FLOAT *y, FLOAT
				*x) 
{
  parms_arms_data data;
  ilutptr ilus;
  int *ipar, schur_start, n;
 
  data = (parms_arms_data)self->data;
  schur_start = data->schur_start;
  ilus   = data->ilus;
  ipar   = data->ipar;

  n = ilus->n;
  if (y == NULL || x == NULL) {
    return 0;
  }
  if (ipar[0] == 0)  {
    invsp(schur_start, ilus, y, x);
  }
  else {
    PARMS_MEMCPY(x, y, n);
    SchLsol(ilus, x);
    SchUsol(ilus, x);

  }

  return 0;
}

static int parms_arms_ascend_vcsr(parms_Operator self, FLOAT *y, FLOAT
				  *x)  
{
  parms_arms_data data;
  p4ptr levmat=NULL, last=NULL;
  ilutptr ilus;
  int *ipar, schur_start, n, nloc, lenB, first;

  data = (parms_arms_data)self->data;
  levmat = data->levmat;
  ilus   = data->ilus;
  schur_start = data->schur_start;
  ipar   = data->ipar;
  n      = data->n;

  if (ipar[0] == 0) {
    Usolp(schur_start, ilus->U, y, x);
    return 0;
  }

  while (levmat) {
    last = levmat;
    levmat = levmat->next;
  }
  levmat = last;

  nloc=levmat->n; 
  lenB=levmat->nB; 
  first = n - nloc; 
  /*-------------------- last level                                 */
  first += lenB; 
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
}

static int parms_arms_getssize_vcsr(parms_Operator self)
{
  parms_arms_data data;

  data = (parms_arms_data)self->data;
  return data->schur_start;
}

static struct parms_Operator_ops parms_arms_sol_vptr = {
  parms_arms_sol_vcsr,
  parms_arms_lsol_vcsr,
  parms_arms_invs_vcsr,
  parms_arms_ascend_vcsr,
  parms_arms_getssize_vcsr,
  parms_arms_nnz,
  arms_free_vcsr,
  arms_view_vcsr
};


int parms_arms_vcsr(parms_Mat self, parms_FactParam param, void *mat,
		    parms_Operator *op)
{
/*---------------------------------------------------------------------
| MULTI-LEVEL BLOCK ILUT PRECONDITIONER.
| ealier  version:  June 23, 1999  BJS -- 
| version2: Dec. 07th, 2000, YS  [reorganized ]
| version 3 (latest) Aug. 2005.  [reorganized + includes ddpq]
+---------------------------------------------------------------------- 
| ON ENTRY:
| ========= 
| ( Amat ) = original matrix A stored in C-style Compressed Sparse
|            Row (CSR) format -- 
|            see LIB/heads.h for the formal definition of this format.
|
| ipar[0:17]  = integer array to store parameters for 
|       arms construction (arms2) 
|
|       ipar[0]:=nlev.  number of levels (reduction processes). 
|                       see also "on return" below. 
| 
|       ipar[1]:= level-reordering option to be used.  
|                 if ipar[1]==0 ARMS uses a block independent set ordering
|                  with a sort of reverse cutill Mc Kee ordering to build 
|                  the blocks. This yields a symmetric ordering. 
|                  in this case the reordering routine called is indsetC
|                 if ipar[1] == 1, then a nonsymmetric ordering is used.
|                  In this case, the B matrix is constructed to be as
|                  diagonally dominant as possible and as sparse as possble.
|                  in this case the reordering routine called is ddPQ.
|                 
|       ipar[2]:=bsize. for indset  Dimension of the blocks. 
|                  bsize is only a target block size. The size of 
|                  each block can vary and is >= bsize. 
|                  for ddPQ: this is only the smallest size of the 
|                  last level. arms will stop when either the number 
|                  of levels reaches nlev (ipar[0]) or the size of the
|                  next level (C block) is less than bsize.
|
|       ipar[3]:=iout   if (iout > 0) statistics on the run are 
|                       printed to FILE *ft
|
|	ipar[4]:= Krylov subspace dimension for last level 
|		    ipar[4] == 0 means only backward/forward solve
|		    is performed on last level.                  
|	ipar[5]:=  maximum # iterations on last level
|
|       ipar[6-9] NOT used [reserved for later use] - set to zero.
| 
| The following set method options for arms2. Their default values can
| all be set to zero if desired. 
|
|       ipar[10-13] == meth[0:3] = method flags for interlevel blocks
|       ipar[14-17] == meth[0:3] = method flags for last level block - 
|       with the following meaning in both cases:
|            meth[0] nonsummetric permutations of  1: yes. affects rperm
|                    USED FOR LAST SCHUR COMPLEMENT 
|            meth[1] permutations of columns 0:no 1: yes. So far this is
|                    USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. 
|                    (so ipar[11] does no matter - enter zero). If 
|                    ipar[15] is one then ILUTP will be used instead 
|                    of ILUT. Permutation data stored in: perm2. 
|            meth[2] diag. row scaling. 0:no 1:yes. Data: D1
|            meth[3] diag. column scaling. 0:no 1:yes. Data: D2
|       all transformations related to parametres in meth[*] (permutation, 
|       scaling,..) are applied to the matrix before processing it 
| 
| ft       =  file for printing statistics on run
|
| droptol  = Threshold parameters for dropping elements in ILU 
|            factorization.
|            droptol[0:4] = related to the multilevel  block factorization
|            droptol[5:6] = related to ILU factorization of last block.
|            This flexibility is more than is really needed. one can use
|            a single parameter for all. it is preferable to use one value
|            for droptol[0:4] and another (smaller) for droptol[5:6]
|            droptol[0] = threshold for dropping  in L [B]. See piluNEW.c:
|            droptol[1] = threshold for dropping  in U [B].
|            droptol[2] = threshold for dropping  in L^{-1} F 
|            droptol[3] = threshold for dropping  in E U^{-1} 
|            droptol[4] = threshold for dropping  in Schur complement
|            droptol[5] = threshold for dropping  in L in last block
|              [see ilutpC.c]
|            droptol[6] = threshold for dropping  in U in last block
|              [see ilutpC.c]
|             This provides a rich selection - though in practice only 4
|             parameters are needed [which can be set to be the same 
              actually] -- indeed it makes sense to take
|             droptol[0] = droptol[1],  droptol[2] = droptol[3], 
|             and droptol[4] = droptol[5]
|
| lfil     = lfil[0:6] is an array containing the fill-in parameters.
|            similar explanations as above, namely:
|            lfil[0] = amount of fill-in kept  in L [B]. 
|            lfil[1] = amount of fill-in kept  in U [B].
|            lfil[2] = amount of fill-in kept  in E L\inv 
|            lfil[3] = amount of fill-in kept  in U \inv F
|            lfil[4] = amount of fill-in kept  in S    .
|            lfil[5] = amount of fill-in kept  in L_S  .
|            lfil[6] = amount of fill-in kept  in U_S 
|             
| tolind   = tolerance parameter used by the indset function. 
|            a row is not accepted into the independent set if 
|            the *relative* diagonal tolerance is below tolind.
|            see indset function for details. Good values are 
|            between 0.05 and 0.5 -- larger values tend to be better
|            for harder problems.
| 
| ON RETURN:
|=============
|
| (PreMat)  = arms data structure which consists of two parts:
|             levmat and ilsch. 
|
| ++(levmat)= permuted and sorted matrices for each level of the block 
|             factorization stored in PerMat4 struct. Specifically
|             each level in levmat contains the 4 matrices in:
|
|
|            |\         |       |
|            |  \   U   |       |
|            |    \     |   F   |
|            |  L   \   |       |
|            |        \ |       |
|            |----------+-------|
|            |          |       |
|            |    E     |       |
|            |          |       |
|            
|            plus a few other things. See LIB/heads.h for details.
|
| ++(ilsch) = IluSpar struct. If the block of the last level is:
|
|                        |  B    F |
|                  A_l = |         | 
|                        |  E    C |
|
|             then IluSpar contains the block C and an ILU
|             factorization (matrices L and U) for the last 
|             Schur complement [S ~ C - E inv(B) F ]
|             (modulo dropping) see LIB/heads.h for details.
|
| ipar[0]   = number of levels found (may differ from input value) 
|
+---------------------------------------------------------------------*/
  /* local matrix object */
  parms_Operator newOpt;
  parms_arms_data data;
  parms_vcsr      Amat;
  parms_Map       is;
/*-------------------- function  prototyping  done in LIB/protos.h    */
/*-------------------- move above to protos.h */
  p4ptr levp=NULL, levc=NULL, levn=NULL, levmat=NULL;
  csptr schur=NULL, B=NULL, F=NULL, E=NULL, C=NULL; 
  ilutptr ilsch=NULL; 
/*-------------------- local variables  (initialized)   */
  double *dd1, *dd2;
  double *droptol, tolind;
  int nlev, bsize, iout, ierr = 0;
  int *ipar, *lfil;
  int methL[4], methS[4];
/*--------------------  local variables  (not initialized)   */
  int nA, nB, nC, j, n, ilev, symperm, schur_start, nbnd, i;
  FILE *ft;
/*--------------------    work arrays:    */
  int *iwork, *uwork; 
/*   timer arrays:  */ 
/*   double *symtime, *unstime, *factime, *tottime;*/
/*---------------------------BEGIN ARMS-------------------------------*/
/*   schur matrix starts being original A */ 

/*-------------------- begin                                         */
  n = param->n;						
  nbnd = schur_start = param->schur_start;
  is = self->is;
  if (schur_start == -1) {
    nbnd = is->schur_start;
  }

  if (!param->isalloc) {
    parms_OperatorCreate(&newOpt);
    PARMS_MEMCPY(newOpt->ops, &parms_arms_sol_vptr, 1);
    PARMS_NEW(data);
    PARMS_NEW(data->levmat);
    levmat = data->levmat;
    PARMS_NEW(data->ilus);
    ilsch =  data->ilus;
    newOpt->data = data;
    *op = newOpt;
  }
  else {
    data = (*op)->data;
  }

  Amat = (parms_vcsr)mat;

  /* compute the number of nonzero entries in mat */
  data->nnz_mat = 0;
  for (i = 0; i < Amat->n; i++) {
    data->nnz_mat += Amat->nnzrow[i];
  }

  /* retrieve data from input objects */
  ipar = param->ipar;
  lfil = param->lfil;
  droptol = param->droptol;
  tolind  = param->tolind;
   
  nlev = ipar[0];
  bsize = ipar[2];
  iout = ipar[3];

  if (iout > 0 ) {
    ft = stdout;
  }
  else{
    ft = NULL;
  }

  ierr = 0;
/*---------------------------------------------------------------------
| The matrix (a,ja,ia) plays role of Schur compl. from the 0th level.
+--------------------------------------------------------------------*/
  
  nC = nA = n = Amat->n;
  if (bsize >= n) bsize = n-1;
  levmat->n = n; levmat->nB = 0;
  schur = Amat;
  levc = levmat;
  /*--------------------------------------- */ 
  levc->prev = levc->next = levp = NULL; 
  levc->n = 0;

  memcpy(methL, &ipar[10], 4*sizeof(int));
  memcpy(methS, &ipar[14], 4*sizeof(int));
  /*---------------------------------------------------------------------
    | The preconditioner construction is divided into two parts:
    |   1st part: construct and store multi-level L and U factors;
    |   2nd part: construct the ILUT factorization for the coarsest level
    +--------------------------------------------------------------------*/
  if ( (iout > 0)  && (nlev > 0) ) {
    fprintf(ft,"  \n");
    fprintf(ft,"Level   Total Unknowns    B-block   Coarse set\n");
    fprintf(ft,"=====   ==============    =======   ==========\n");
  }
  /*---------------------------------------------------------------------
    | main loop to construct multi-level LU preconditioner. Loop is on the
    | level ilev. nA is the dimension of matrix A_l at each level.
    +--------------------------------------------------------------------*/
  for (ilev = 0; ilev < nlev; ilev++) {
    /*-------------------- new nA is old nC -- from previous level */
    nA = nC;
    if ( nA <= bsize )  goto label1000;  
    /*-------------------- allocate work space                        */ 
    iwork = (int *) Malloc(nA*sizeof(int), "arms2:2.5" );
    symperm = 0;    /* 0nly needed in cleanP4 */
    if (ipar[1] == 1) 
      uwork = (int *) Malloc(nA*sizeof(int), "arms2:2.6" );
    else{
      symperm = 1;    
      uwork = iwork; 
    }
    /*-------------------- SCALING*/
    dd1 = NULL;
    dd2 = NULL;
    if (methL[2]) {
      dd1 = (double *) Malloc(nA*sizeof(double), "arms2:3" );
      j=roscalC(schur, dd1,1);
      if (j) printf("ERROR in roscalC -  row %d  is a zero row\n",j);
    }

    if (methL[3]) {
      dd2 = (double *) Malloc(nA*sizeof(double), "arms2:4" );
      j=coscalC(schur, dd2,1); 
      if (j) printf("ERROR in coscalC - column %d is a zero column\n",j);
    }
    /*--------------------independent-sets-permutation-------------------
      |  do reordering -- The matrix and its transpose are used.
      +--------------------------------------------------------------------*/
    /* if (SHIFTTOL > 0.0) shiftsD(schur,SHIFTTOL);    */
    if (ipar[1] == 1) 
      PQperm(schur, uwork, iwork, &nB, tolind, nbnd) ; 
    else 
      indsetC (schur, bsize, iwork, &nB, tolind, nbnd) ; 
    /*---------------------------------------------------------------------
      | nB is the total number of nodes in the independent set.
      | nC : nA - nB = the size of the reduced system.
      +--------------------------------------------------------------------*/
    nC = nA - nB;
    nbnd -= nB;
    /*   if the size of B or C is zero , exit the main loop  */
    /*   printf ("  nB %d nC %d \n",nB, nC); */
    if ( nB == 0 || nC == 0 )  goto label1000; 
    /*---------------------------------------------------------------------
      | The matrix for the current level is in (schur).
      | The permutations arrays are in iwork and uwork (row).
      | The routines rpermC, cpermC permute the matrix in place.
      *-----------------------------------------------------------------------*/
    /*   DEBUG : SHOULD THIS GO BEFORE GOTO LABEL1000 ?? */
    rpermC(schur,uwork); 
    cpermC(schur,iwork);
    /*   prtC(schur, ilev) ;   print matrix - debugging */
    /*-----------------------------------------------------------------------
      | If this is the first level, the permuted matrix is stored in 
      | (levc) = (levmat).  Otherwise, the next level is created in (levc).
      +--------------------------------------------------------------------*/
    if (ilev > 0) {
      /*-   delete C matrix of any level except last one (no longer needed) */
      cleanCS(C); 
      levn = (p4ptr) Malloc(sizeof(Per4Mat), "arms2:6" );
      /* levc->prev = levp; */
      levc->next = levn;
      levp = levc;
      levc = levn;
      levc->prev = levp; 
    }
    /*-------------------- p4ptr struct for current schur complement */
    B = (csptr) Malloc(sizeof(SparMat), "arms2:7" );
    E = (csptr) Malloc(sizeof(SparMat), "arms2:8" );
    F = (csptr) Malloc(sizeof(SparMat), "arms2:9" );
    C = (csptr) Malloc(sizeof(SparMat), "arms2:10" );
    csSplit4(schur, nB, nC, B, F, E, C);
    setupP4(levc, nB, nC, F, E);
    /*--------------------     copy a few pointers       ---- */      
    levc->perm  = iwork;
    levc->rperm = uwork; 
    levc->symperm = symperm;
    levc->D1=dd1;
    levc->D2=dd2; 
    /*---------------------------------------------------------------------
      | a copy of the matrix (schur) has been permuted. Now perform the 
      | block factorization: 
      |
      | | B   F |       | L       0 |     | U  L^-1 F |
      | |       |   =   |           |  X  |           | = L x U
      | | E   C |       | E U^-1  I |     | 0    A1   |
      |   
      | The factors E U^-1 and L^-1 F are discarded after the factorization.
      |
      +--------------------------------------------------------------------*/ 
    if (iout > 0)
      fprintf(ft,"%3d %13d %13d %10d\n", ilev+1,nA,nB,nC);
    /*---------------------------------------------------------------------
      | PILUT constructs one level of the block ILU fact.  The permuted matrix
      | is in (levc).  The L and U factors will be stored in the p4mat struct.
      | destroy current Schur  complement - no longer needed  - and set-up new
      | one for next level...
      +--------------------------------------------------------------------*/
    cleanCS(schur);
    schur = (csptr) Malloc(sizeof(SparMat), "arms2:11" ); 
    setupCS(schur, nC,1);
    /*----------------------------------------------------------------------
      | calling PILU to construct this level block factorization
      | ! core dump in extreme case of empty matrices.
      +----------------------------------------------------------------------*/
    ierr = pilu(levc, B, C, droptol, lfil, schur) ;
    /* prtC(levc->L, ilev) ; */
    if (ierr) { 
      fprintf(ft," ERROR IN  PILU  -- IERR = %d\n", ierr);
      return(1);
    }
    cleanCS(B); 
  }
  /*---------------------------------------------------------------------
    |   done with the reduction. Record the number of levels in ipar[0] 
    |**********************************************************************
    +--------------------------------------------------------------------*/
 label1000:
  /* printf (" nnz_Schur %d \n",cs_nnz (schur)); */
  levc->next = NULL;
  ipar[0] = ilev;
  data->nlev = ilev;  
  data->n = n; 
  nC = schur->n;
  setupILUT(ilsch,nC); 
  /*--------------------------------------------------------------------*/
  /* define C-matrix (member of ilsch) to be last C matrix */ 
  if (ilev > 0) ilsch->C=C; 
  /*-------------------- for ilut fact of schur matrix */
  /*  SCALING  */

  ilsch->D1 = NULL;
  if (methS[2]) {
    ilsch->D1 = (double *) Malloc(nC*sizeof(double), "arms2:iluschD1" );
    j=roscalC(schur, ilsch->D1, 1); 
    if (j) printf("ERROR in roscalC - row %d is a zero row\n",j);
  }

  ilsch->D2  = NULL;
  if (methS[3]) {
    ilsch->D2 = (double *) Malloc(nC*sizeof(double), "arms2:iluschD1" );
    j =coscalC(schur, ilsch->D2, 1);  
    if (j) printf("ERROR in coscalC - column %d is a zero column\n",j);
  }

  /*---------------------------------------------------------------------
    |     get ILUT factorization for the last reduced system.
    +--------------------------------------------------------------------*/
  uwork = NULL;
  iwork = NULL;
  if (methS[0]) { 
    iwork = (int *) Malloc(nC*sizeof(int), "arms2:3" );
    uwork = (int *) Malloc(nC*sizeof(int), "arms2:3.5" );
    tolind = 0.0; 
    PQperm(schur, uwork, iwork, &nB, tolind, nbnd) ; 
    rpermC(schur,uwork); 
    cpermC(schur,iwork);
  }
  ilsch->rperm = uwork; 
  ilsch->perm  = iwork;

  /*   printf("  lf : %d  %d  %d  %d  %d  %d  %d  \n",lfil[0],  
       lfil[1], lfil[2], lfil[3], lfil[4], lfil[5], lfil[6]) ; */
   
  ilsch->perm2 = NULL; 

  if (methS[1] == 0)
    ierr = ilutD(schur, droptol, lfil, ilsch);
  else {
    ilsch->perm2 = (int *) Malloc(nC*sizeof(int), "arms2:ilutpC" );
    for (j=0; j<nC; j++)
      ilsch->perm2[j] = j;
    ierr = ilutpC(schur, droptol, lfil, PERMTOL, nC, ilsch);
  }
  /*---------- OPTIMIZATRION: NEED TO COMPOUND THE TWO
    RIGHT PERMUTATIONS -- CHANGES HERE AND IN 
    USCHUR SOLVE ==  compound permutations */     
  if (ierr) {
    fprintf(ft," ERROR IN  ILUT -- IERR = %d\n", ierr); 
    return(1); 
  }
  /* Last Schur complement no longer needed */
  cleanCS(schur);
  data->nnz_prec = nnz_arms(data);
  data->ind = n - ilsch->n;
  if (ilev) {
    data->schur_start = n - ilsch->n;
  }
  else {
    is = self->is;
    data->schur_start = is->schur_start;
  }

  PARMS_MEMCPY(data->ipar,   param->ipar, 18);
  PARMS_MEMCPY(data->pgfpar, param->pgfpar, 2);
  return 0;
}/*-----end-of-ARMS2----------------------------------------------------
   +--------------------------------------------------------------------*/

