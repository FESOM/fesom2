#include <math.h>
#include "include/parms_mat_impl.h"
#include "include/parms_opt_impl.h"
#include "DDPQ/protos.h"

#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

typedef struct parms_ilu_data {
  parms_vcsr L;
  parms_vcsr U;
  int schur_start;
  int n;
  int *nnzschur;
  int nnz_mat;
  int nnz_prec;
} *parms_ilu_data;

/*
int qsplitC(FLOAT *a, int *ind, int n, int ncut);
void qqsort(int *ja, FLOAT *ma, int left, int right);
*/

static int parms_ilu_free(parms_Operator *self)
{
  parms_ilu_data data;
  int n, nnz, i, *pj;
  FLOAT *pa;

  data = (parms_ilu_data)(*self)->data;
  n = data->L->n;

  for (i = 0; i < n; i++) {
    nnz = data->L->nnzrow[i];
    pj  = data->L->pj[i];
    pa  = data->L->pa[i];
    if (nnz) {
      PARMS_FREE(pj);
      PARMS_FREE(pa);
    }
  }
  PARMS_FREE(data->L->nnzrow);
  PARMS_FREE(data->L->pj);
  PARMS_FREE(data->L->pa);
  PARMS_FREE(data->L);

  n = data->U->n;
  for (i = 0; i < n; i++) {
    nnz = data->U->nnzrow[i];
    pj  = data->U->pj[i];
    pa  = data->U->pa[i];
    if (nnz) {
      PARMS_FREE(pj);
      PARMS_FREE(pa);
    }
  }

  PARMS_FREE(data->U->nnzrow);
  PARMS_FREE(data->U->pj);
  PARMS_FREE(data->U->pa);
  PARMS_FREE(data->U);

  if (data->n - data->schur_start){
    PARMS_FREE(data->nnzschur);
  }

  PARMS_FREE(data);
  return 0;
}
    
static int parms_ilu_view(parms_Operator self, parms_Viewer v)
{
  parms_ilu_data data;
  int i, j, n, nnz, *pj;
  FLOAT *pa;
  FILE *fp;

  parms_ViewerGetFP(v, &fp);
  data = (parms_ilu_data)self->data;
  n = data->L->n;
  fprintf(fp, "L part of the matrix:\n");
  fprintf(fp, "n = %d\n", n);
#if defined(DBL_CMPLX)  
  for (i = 0; i < n; i++) {
    nnz = data->L->nnzrow[i];
    pj  = data->L->pj[i];
    pa  = data->L->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);    
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f  %f) ", i, pj[j], creal(pa[j]), cimag(pa[j]));
    }
  }
  n = data->U->n;
  fprintf(fp, "U part of the matrix:\n");
  fprintf(fp, "n = %d\n", n);
  for (i = 0; i < n; i++) {
    nnz = data->U->nnzrow[i];
    pj  = data->U->pj[i];
    pa  = data->U->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f  %f) ", i, pj[j], creal(pa[j]), cimag(pa[j]));
    }
  }
#else
  for (i = 0; i < n; i++) {
    nnz = data->L->nnzrow[i];
    pj  = data->L->pj[i];
    pa  = data->L->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);    
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f) ", i, pj[j], pa[j]);
    }
  }
  n = data->U->n;
  fprintf(fp, "U part of the matrix:\n");
  fprintf(fp, "n = %d\n", n);
  for (i = 0; i < n; i++) {
    nnz = data->U->nnzrow[i];
    pj  = data->U->pj[i];
    pa  = data->U->pa[i];
    fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
    for (j = 0; j < nnz; j++) {
      fprintf(fp, "(%d,%d,%f) ", i, pj[j], pa[j]);
    }
  }
#endif
  parms_ViewerStoreFP(v, fp);
  return 0;
}

static void parms_ilu_nnz_vcsr(parms_Operator self, int *nnz_mat, int *nnz_pc)
{
  parms_ilu_data data;

  data = (parms_ilu_data)self->data;
  *nnz_mat = data->nnz_mat;
  *nnz_pc  = data->nnz_prec;
}

static int parms_ilu_lsol_vcsr(parms_Operator self, FLOAT *y, FLOAT
			       *x) 
{
  parms_ilu_data data;
  parms_vcsr L;
  int n, i, j, nnz, *pj, start;
  FLOAT *pa, t;

  data = (parms_ilu_data)self->data;
  start = data->schur_start;
  L = data->L;
  n = L->n;
  for (i = 0; i < start; i++) {
    nnz = L->nnzrow[i];
    pj  = L->pj[i];
    pa  = L->pa[i];
    t = y[i];
    for (j = 0; j < nnz; j++) {
      t -= pa[j] * x[pj[j]];
    }
    x[i] = t;
  }

  for (i = start; i < n; i++) {
    pj  = L->pj[i];
    pa  = L->pa[i];
    t   = y[i];
    nnz = data->nnzschur[i-start];
    for (j = 0; j < nnz; j++) {
      t -= pa[j] * x[pj[j]];
    }
    x[i] = t;
  }

  return 0;
}

static int parms_ilu_invs_vcsr(parms_Operator self,  FLOAT *y, FLOAT
			       *x) 
{
  parms_ilu_data data;
  parms_vcsr L, U;
  int n, i, j, nnz, start, *pj, nnzS;
  FLOAT *pa, t, diag;

  data = (parms_ilu_data)self->data;
  start = data->schur_start;
  L = data->L;
  n = L->n;
  for (i = start; i < n; i++) {
    nnz  = L->nnzrow[i];
    pj   = L->pj[i];
    pa   = L->pa[i];
    t    = y[i-start];
    nnzS = data->nnzschur[i-start];
    for (j = nnzS; j < nnz; j++) {
      t -= pa[j] * x[pj[j]-start];
    }
    x[i-start] = t;
  }

  U = data->U;
  for (i = n-1; i >= start; i--) {
    nnz = U->nnzrow[i];
    pj  = U->pj[i];
    pa  = U->pa[i];
    t   = x[i-start];
    diag = pa[0];
    for (j = 1; j < nnz; j++) {
      t -= pa[j] * x[pj[j]-start];
    }
    x[i-start] = t*diag;
  }
  return 0;
}

static int parms_ilu_getu_vcsr(parms_Operator self, void **mat)
{
  parms_ilu_data data;
  parms_vcsr U, S;
  int i, n, start, nnz;
  
  data = (parms_ilu_data)self->data;
  start = data->schur_start;
  U = data->U;
  *mat = U;
  return 0;
}


static int parms_ilu_ascend_vcsr(parms_Operator self, FLOAT *y, FLOAT
				 *x)  
{
  parms_ilu_data data;
  parms_vcsr U;
  int i, j, nnz, start, *pj;
  FLOAT *pa, t, diag;

  data = (parms_ilu_data)self->data;

  U = data->U;
  start = data->schur_start;

  for (i = start-1; i >= 0; i--) {
    nnz = U->nnzrow[i];
    pj  = U->pj[i];
    pa  = U->pa[i];
    t = y[i];
    diag = pa[0];
    for (j = 1; j < nnz; j++) {
      t -= pa[j] * x[pj[j]];
    }
    x[i] = t*diag;
  }
  return 0;
}

static int parms_ilu_sol_vcsr(parms_Operator self, FLOAT *y,
			      FLOAT *x)  
{
  FLOAT *pa;
  parms_ilu_data data;
  parms_vcsr L, U;
  int i, j, n, *pj, nnz;

  data = (parms_ilu_data)self->data;
  L = data->L;
  U = data->U;
  n = L->n;

  for (i = 0; i < n; i++) {
    // nnz = L->nnzrow[i];
    pj  = L->pj[i];
    pa  = L->pa[i];
    x[i] = y[i];

    for (j = 0; j < L->nnzrow[i]; j++) {
      x[i] -= pa[j] * x[pj[j]];
    }
  }

  for (i = n-1; i >= 0; i--) {
    //    nnz = U->nnzrow[i];
    pj  = U->pj[i];
    pa  = U->pa[i];

    for (j = 1; j < U->nnzrow[i]; j++) {
      x[i] -= pa[j] * x[pj[j]];
    }
    x[i] *= pa[0];  // pa[0] = diag
  }

  return 0;
}

static int parms_ilu_getssize_vcsr(parms_Operator self)
{
  parms_ilu_data data;

  data = (parms_ilu_data)self->data;
  return data->schur_start;
}
  
static struct parms_Operator_ops parms_ilu_sol_vptr = {
  parms_ilu_sol_vcsr,
  parms_ilu_lsol_vcsr,
  parms_ilu_invs_vcsr,
  parms_ilu_getu_vcsr,
  parms_ilu_ascend_vcsr,
  parms_ilu_getssize_vcsr,
  parms_ilu_nnz_vcsr,
  parms_ilu_free,
  parms_ilu_view
};

int parms_ilu0_vcsr(parms_Mat self, parms_FactParam param, void *mat,
		    parms_Operator *op) 
{
  /*
   * it is assumed that the the elements in the input matrix are
   * stored in such a way that in each row the lower part comes first
   * and then the upper part. To get the correct ILU factorization, it
   * is also necessary to have the elements of L sorted by increasing
   * column number. It may therefore be necessary to sort the elements
   * of a, ja, ia prior to calling ilu0. This can be achieved by
   * transposing the matrix twice using csrcsc.
   *------------------------------------------------------------------*/
  parms_Operator newOpt;
  parms_ilu_data data;
  parms_vcsr     mat_vcsr;
  parms_Map      is;
  int n, start, schur_start, m, *iw, ierr, ii, jj, i, j, k, nnz, lenl, lenu;
  int *pj, *rowj, *rowjj, nloc, jcol, jpos, jrow, jw;
  FLOAT *pa, *rowm, *rowmm, t1, shift, rnrm;
//  int flag;
  BOOL found;
  int rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  is = self->is;
  n = param->n;
  jw = n; /* initialize jw */
  nloc = parms_MapGetLocalSize(is);
  start = param->start;					
  schur_start = param->schur_start;
  if (schur_start == -1) {
    schur_start = n;
  }

  if (!param->isalloc) {
    parms_OperatorCreate(&newOpt);
    PARMS_MEMCPY(newOpt->ops, &parms_ilu_sol_vptr, 1);		
    PARMS_NEW(data);						
    PARMS_NEW(data->L);					
    data->L->n = n;						
    PARMS_NEWARRAY(data->L->nnzrow, n);			
    PARMS_NEWARRAY(data->L->pj, n);				
    PARMS_NEWARRAY(data->L->pa, n);				
    PARMS_NEW(data->U);					
    data->U->n = n;						
    PARMS_NEWARRAY(data->U->nnzrow, n);			
    PARMS_NEWARRAY(data->U->pj,     n);			
    PARMS_NEWARRAY(data->U->pa,     n);
    data->nnz_mat  = 0;
    data->nnz_prec = 0;
    param->isalloc = true;
    
    newOpt->data = data;
    *op = newOpt;
  }
  else {
    data = (*op)->data;
  }

  mat_vcsr = (parms_vcsr)mat;

  PARMS_NEWARRAY(iw, n);
  for(i = 0; i < n; i++) {
    iw[i] = -1;
  }
  /* main loop */
  ierr = 0;
  m = mat_vcsr->n;

  /* sort the column indices of the entries received from neighbouring
   * processors in ascending order */
  for (ii = nloc; ii < m+start; ii++) {
    nnz = mat_vcsr->nnzrow[ii];
    pj  = mat_vcsr->pj[ii];
    pa  = mat_vcsr->pa[ii];
    qqsort(pj, pa, 0, nnz-1);
  }

  for(ii = 0; ii < m; ii++) {
    nnz = mat_vcsr->nnzrow[ii];
    pj  = mat_vcsr->pj[ii];
    pa  = mat_vcsr->pa[ii];

    /* calculating the number of entries in L and U */
    lenl = 0;
    lenu = 0;
    if (ii+start < schur_start) {
      for(j = 0; j < nnz; j++) {
	jcol = pj[j];
	if(jcol < ii+start) {
	  ++lenl;
	}
	else {
	  ++lenu;
	}
      }
    }
    else {
      for(j = 0; j < nnz; j++) {
	jcol = pj[j];
	if(jcol < schur_start) {
	  ++lenl;
	}
	else {
	  ++lenu;
	}
      }
    }

    data->L->nnzrow[ii+start] = lenl;
      
    if(lenl > 0) {
      PARMS_NEWARRAY(data->L->pj[ii+start], lenl);
      PARMS_NEWARRAY(data->L->pa[ii+start], lenl);
    }
    if(lenu == 0)
    {
      if(lenl == 0){
         printf("myid %d: Zero row encountered at row %d \n", rank, ii);
         MPI_Abort(MPI_COMM_WORLD, 10);
      }
    }
    
    data->U->nnzrow[ii+start] = lenu;
/* Assume diagonal element */
    PARMS_NEWARRAY(data->U->pj[ii+start], lenu+1);
    PARMS_NEWARRAY(data->U->pa[ii+start], lenu+1);
    
    lenl = 0;
    lenu = 0;
    shift = 0.0;
    rnrm = 0.0;
    /* copy row ii into L and U part */
    if (ii+start < schur_start) {
      lenu = 1; // Assume diagonal element   
      found = false;    
      data->U->pa[ii+start][0] = 0;	  
	data->U->pj[ii+start][0] = ii+start;
	iw[ii+start] = ii + start;      
                 
      for(j = 0; j < nnz; j++) {
	jcol = pj[j];
	rnrm += ABS_VALUE(pa[j]);
	shift = ABS_VALUE(shift) > ABS_VALUE(pa[j]) ? shift : pa[j]; 
	if(jcol < ii+start) {
	  data->L->pj[ii+start][lenl] = jcol;
	  data->L->pa[ii+start][lenl] = pa[j];
	  iw[jcol] = lenl;
	  ++lenl;
	}
	else if(jcol == ii + start){
	/* Diagonal element */
	  data->U->pa[ii+start][0] = pa[j];	  
	  found = true;
	}	  	  
	else {
	  data->U->pj[ii+start][lenu] = jcol;
	  data->U->pa[ii+start][lenu] = pa[j];
	  iw[jcol] = ii + start+ lenu;
	  ++lenu;
	 }
      }
/* create room for diagonal element, if none.
 * This just ensures that the fill-factor is 1.
 */      
      shift /= rnrm;
      if(found == false){
        data->U->nnzrow[ii+start] += 1;
	}
/* perturb zero or small diagonals - this is done later*/      	
//	else if(ABS_VALUE(data->U->pa[ii+start][0]) <= DBL_EPSILON)
//	      data->U->pa[ii+start][0] = 0.0;//shift;    
    }
    else {
      for(j = 0; j < nnz; j++) {
	jcol = pj[j];
	if(jcol < schur_start) {
	  data->L->pj[ii+start][lenl] = jcol;
	  data->L->pa[ii+start][lenl] = pa[j];
	  iw[jcol] = lenl;
	  ++lenl;
	}
	else {
	  data->U->pj[ii+start][lenu] = jcol;
	  data->U->pa[ii+start][lenu] = pa[j];
	  iw[jcol] = schur_start + lenu;
	  ++lenu;
	}
      }
    }

    /* exit if diagonal element is reached. */
    rowj = data->L->pj[ii+start];
    rowm = data->L->pa[ii+start];
    for(j = 0; j < lenl; j++) {
/*----------------------------------------------------------------------------
 *  in order to do the elimination in the correct order we must select the
 *  smallest column index among rowj[k], k = j+1, ..., lenl
 *--------------------------------------------------------------------------*/
      jrow = rowj[j];
      jpos = j;
/*---------- determine smallest column index */
      for( k = j + 1; k < lenl; k++ ) {
        if( rowj[k] < jrow ) {
          jrow = rowj[k];
          jpos = k;
        }
      }
      if( jpos != j ) {
/*---------- swaps */
        jcol = rowj[j];
        rowj[j] = rowj[jpos];
        rowj[jpos] = jcol;
        iw[jrow] = j;
        iw[jcol]  = jpos;
        t1 = rowm[j];
        rowm[j] = rowm[jpos];
        rowm[jpos] = t1;
      }    
    
 /* Eliminate previous rows: jrow now contains the smallest 
    column index in L part
 */   
//      jrow = rowj[j];
      rowjj = data->U->pj[jrow];
      rowmm = data->U->pa[jrow];
      nnz = data->U->nnzrow[jrow];
      rowm[j] *= rowmm[0];
      t1 = rowm[j];
      /* perform linear combination */
      if (ii+start < schur_start) {
	for(jj = 1; jj < nnz; jj++) {
	  jw = iw[rowjj[jj]];
	  if(jw != -1) {
	    if(jw < ii+start) {
	      rowm[jw] -= t1*rowmm[jj];
	    }
	    else {
	      data->U->pa[ii+start][jw-ii-start] -= t1*rowmm[jj];
	    }
	  }
	}
      }
      else {
	for (jj = 1; jj < nnz; jj++) {
	  if(jw < schur_start) {
	    rowm[jw] -= t1*rowmm[jj];
	  }
	  else {
	    data->U->pa[ii+start][jw-schur_start] -= t1*rowmm[jj];
	  }
	}
      }
    }

    if (ii+start < schur_start) {
        t1 = data->U->pa[ii+start][0];
/* Check for zero diagonal */    
        if(ABS_VALUE(t1) <= DBL_EPSILON) {
/*          ierr = ii+1; 
            printf("myid %d: Zero diagonal encountered at row %d \n", rank, ii);
            MPI_Abort(MPI_COMM_WORLD, 10);
*/
/* Perturb zero diagonal to be (maximum elimination coeff) */
            t1 = shift;      
         }
      /* invert and store diagonal element. */
         data->U->pa[ii+start][0] = 1.0/t1;
    }
    /* reset pointer iw to -1 */
    rowj = data->L->pj[ii+start];
    for(i = 0; i < lenl; i++) {
      iw[rowj[i]] = -1;
    }
    rowj = data->U->pj[ii+start];
    for(i = 0; i < lenu; i++) {
      iw[rowj[i]] = -1;
    }
  }

  data->schur_start = is->schur_start;
  data->n = n;

    if (n-data->schur_start) {
      FLOAT *w;
      PARMS_NEWARRAY(w, n);
      PARMS_NEWARRAY(data->nnzschur, n-data->schur_start);
      for (i = data->schur_start; i < n; i++) {
	lenl = data->L->nnzrow[i];
	rowj = data->L->pj[i];
	rowm = data->L->pa[i];
	jw   = 0;

	for (j = 0; j < lenl; j++) {
	  if (rowj[j] < data->schur_start) {
	    ++jw;
	  }
	}
	
	data->nnzschur[i-data->schur_start] = jw;

	ii = 0;
	jj = jw;
	for (j = 0; j < lenl; j++) {
	  if (rowj[j] < data->schur_start) {
	    iw[ii]  = rowj[j];
	    w[ii++] = rowm[j];
	  }
	  else {
	    iw[jj]  = rowj[j];
	    w[jj++] = rowm[j];
	  }
	}
	PARMS_MEMCPY(rowj, iw, lenl);
	PARMS_MEMCPY(rowm, w,  lenl);
      }
      PARMS_FREE(w);
    }
  PARMS_FREE(iw);

  /* compute the number of nonzeros in matrix */
  for (i = 0; i < m; i++) {
    data->nnz_mat += mat_vcsr->nnzrow[i];
  }

  /* compute the number of nonzeros in pc */
  for (i = 0; i < m; i++) {
    data->nnz_prec += data->L->nnzrow[i];
    data->nnz_prec += data->U->nnzrow[i];
  }
  return ierr;
}

#define MIN(a, b) (a) > (b) ? (b) : (a)
#define MAX_ROWS_PER_COL 20 
int parms_iluk_vcsr(parms_Mat self, parms_FactParam param, void *mat,
		    parms_Operator *op) 
{
  /*
   * ON RETURN
   *===========
   *  ierr    =  integer. Error message with the following meaning.
   *             ierr  = 0    --> successful return.
   *             ierr  > 0   --> zero pivot encountered at step
   *                              number ierr.
   *             ierr  = -1   --> Error. input matrix may be wrong.
   *                              (The elimination process has
   *                              generated a row in L or U whose
   *                              length is greater than n)
   *             ierr  = -2   --> Illegal value for lfil.
   *             ierr  = -3   --> zero row encountered in A or U
   *------------------------------------------------------------------*/
  parms_Operator newOpt;
  parms_ilu_data data;
  parms_vcsr     mat_vcsr;
  parms_Map      is;
  int n, start, schur_start, m, nnz, lfil, *jw, **levs, n2, *pj, ierr;
  int lenl, lenu, ii, jj, i, j, k, jcol, jpos, jrow, *rowj, *lev, jlev;
  FLOAT *w, *pa, *rowm, t, s, fact, shift, rnrm;
//  int flag;
//  BOOL found;
  int rank;
  int lenl_all, lenu_all, max_entr;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  is = self->is;
  n = param->n;						
  start = param->start;					
  schur_start = param->schur_start;
  if (schur_start == -1) {
    schur_start = n;
  }

  if (!param->isalloc) {			
// global data arrays, F77 maximum size... 
    max_entr = n*MAX_ROWS_PER_COL;
    parms_OperatorCreate(&newOpt);
    PARMS_MEMCPY(newOpt->ops, &parms_ilu_sol_vptr, 1);		
    PARMS_NEW(data);						
    PARMS_NEW(data->L);					
    data->L->n = n;						
    PARMS_NEWARRAY(data->L->nnzrow, n);			
    PARMS_NEWARRAY(data->L->pj,     n);				
    PARMS_NEWARRAY(data->L->pa,     n);				
// allocate global data arrays 
    PARMS_NEWARRAY(data->L->pj_data, max_entr);				
    PARMS_NEWARRAY(data->L->pa_data, max_entr);					
    PARMS_NEW(data->U);					
    data->U->n = n;						
    PARMS_NEWARRAY(data->U->nnzrow, n);			
    PARMS_NEWARRAY(data->U->pj,     n);			
    PARMS_NEWARRAY(data->U->pa,     n);
// allocate global data arrays 
    PARMS_NEWARRAY(data->U->pj_data, max_entr);				
    PARMS_NEWARRAY(data->U->pa_data, max_entr);		
    data->nnz_mat  = 0;
    data->nnz_prec = 0;
    param->isalloc = true;

    newOpt->data = data;
    *op = newOpt;
  }
  else {
    data = (*op)->data;
  }

  mat_vcsr = (parms_vcsr)mat;

  lfil = param->ipar[0];

  if(lfil < 0) {
    ierr = -2;  
    printf("myid %d: lfil parameter must be > zero: ierr %d \n", rank, ierr);
    MPI_Abort(MPI_COMM_WORLD, 11);
  }

  /* malloc memory for working arrays levs, w and jw */
  PARMS_NEWARRAY(levs, n);
  PARMS_NEWARRAY(w,    n);
  PARMS_NEWARRAY(jw,   3*n);

  /* main loop */
  n2 = n << 1;
  for(i = n; i < n2; i++) {
    jw[i] = -1;
    jw[n+i] = 0;
  }
  m = mat_vcsr->n;
  ierr = 0;
  lenl_all = 0;
  lenu_all = 0;
  for(ii = 0; ii < m; ii++) {
    lenl = 0;
    lenu = 0;
    shift = 0.0;
    rnrm = 0.0;
    nnz = mat_vcsr->nnzrow[ii];
    pj  = mat_vcsr->pj[ii];
    pa  = mat_vcsr->pa[ii];
    if (start+ii < schur_start) {
      lenu = 1;
//      found = false;
      jw[ii+start] = ii + start;
      w[ii+start] = 0.0;
      jw[n+ii+start] = ii + start;

      for(j = 0; j < nnz; j++) {
	jcol = pj[j];
	t    = pa[j];
        if(ABS_VALUE(t) <=  DBL_EPSILON) continue;
      rnrm += ABS_VALUE(t);
	shift = ABS_VALUE(shift) > ABS_VALUE(t) ? shift : t;       
	if(jcol < ii+start) {
	  jw[lenl] = jcol;
	  w[lenl] = t;
	  jw[n2+lenl] = 0;
	  jw[n+jcol] = lenl;
	  ++lenl;
	}
	else if(jcol == ii+start) {
	  w[ii+start] = t;
	  jw[n2+ii+start] = 0;
//	  found = true;
	}
	else {
	  jpos = ii + start + lenu;
	  jw[jpos] = jcol;
	  w[jpos] = t;
	  jw[n2+jpos] = 0;
	  jw[n+jcol] = jpos;
	  ++lenu;
	}
    }
/* perturb zero or small diagonals */      
      shift /= rnrm;
/*      if(found == false){
	  w[ii+start] = shift;	  
	  jw[n2+ii+start] = 0;
	}*/        
    }
    else {
      for(j = 0; j < nnz; j++) {
	jcol = pj[j];
	t    = pa[j];
        if(ABS_VALUE(t) <= DBL_EPSILON) continue;
	if(jcol < schur_start) {
	  jw[lenl] = jcol;
	  w[lenl] = t;
	  jw[n2+lenl] = 0;
	  jw[n+jcol] = lenl;
	  ++lenl;
	}
	else {
	  jpos = schur_start + lenu;
	  jw[jpos] = jcol;
	  w[jpos] = t;
	  jw[n2+jpos] = 0;
	  jw[n+jcol] = jpos;
	  ++lenu;
	}
      }
    }
  
    /* eliminate previous rows */
    for(jj = 0; jj < lenl; jj++) {
      /* in order to do the elimination in the correct order we must
	 select the smallest column index among jw[k], k=jj,...lenl */
      jrow = jw[jj];
      k = jj;
      /* determine smallest column index */
      for(j = jj+1; j < lenl; j++) {
	if(jw[j] < jrow) {
	  jrow = jw[j];
	  k = j;
	}
      }
      if(k != jj) {
	/* exchange in jw */
	j = jw[jj];
	jw[jj] = jw[k];
	jw[k] = j;
	/* exchange in jw(n+ (pointers/ nonzero indicator). */
	jw[n+jrow] = jj;
	jw[n+j] = k;
	/* exchange in jw(n2 + (levels) */
	j = jw[n2+jj];
	jw[n2+jj] = jw[n2+k];
	jw[n2+k] = j;
	/* exchange in w */
	s = w[jj];
	w[jj] = w[k];
	w[k] = s;
      }
      /* zero out element in row by resetting jw[n+jrow] to zero */
      jw[n+jrow] = -1;
      /* get the multiplier for row to be eliminated (jrow) + its
	 level */
      rowj = data->U->pj[jrow];
      rowm = data->U->pa[jrow];
      nnz =  data->U->nnzrow[jrow];
      lev = levs[jrow];
      fact = w[jj]*rowm[0];
      jlev = jw[n2+jj];
      if(jlev > lfil) continue;
      /* combine current row and row jrow */
      if (ii+start < schur_start) {
	for(k = 1; k < nnz; k++) {
	  s = fact*rowm[k];
	  jcol = rowj[k];
	  jpos = jw[n+jcol];
	  if(jcol >= ii+start) {
	    /* dealing with upper part */
	    if(jpos == -1) {
	      /* this is a fill-in element */
	      i = ii + start + lenu;
	      jw[i] = jcol;
	      jw[n+jcol] = i;
	      w[i] = -s;
	      jw[n2+i] = jlev+lev[k]+1;
	      ++lenu;
	      if(lenu > n) {
		  ierr = -1;
              printf("myid %d: Error in ILUK (nnz(U) > n): ierr = %d \n", rank, ierr);
              MPI_Abort(MPI_COMM_WORLD, 12);		
	      }
	    }
	    else {
	      /* this is not a fill-in element */
	      w[jpos] -= s;
	      jw[n2+jpos] = MIN(jw[n2+jpos], jlev+lev[k]+1);
	    }
	  }
	  else {
	    /* dealing with lower part */
	    if(jpos == -1) {
	      /* this is a fill-in element */
	      jw[lenl] = jcol;
	      jw[n+jcol] = lenl;
	      w[lenl] = -s;
	      jw[n2+lenl] = jlev+lev[k]+1;
	      ++lenl;
	      if(lenl > n) {
		ierr = -1;
		break;
	      }
	    }
	    else {
	      /* this is not a fill-in element */
	      w[jpos] -= s;
	      jw[n2+jpos] = MIN(jw[n2+jpos], jlev+lev[k]+1);
	    }
	  }
	}
      }
      else {
	for(k = 1; k < nnz; k++) {
	  s = fact*rowm[k];
	  jcol = rowj[k];
	  jpos = jw[n+jcol];
	  if(jcol >= schur_start) {
	    /* dealing with upper part */
	    if(jpos == -1) {
	      /* this is a fill-in element */
	      i = schur_start + lenu;
	      jw[i] = jcol;
	      jw[n+jcol] = i;
	      w[i] = -s;
	      jw[n2+i] = jlev+lev[k]+1;
	      ++lenu;
	      if(lenu > n) {
		ierr = -1;
		break;
	      }
	    }
	    else {
	      /* this is not a fill-in element */
	      w[jpos] -= s;
	      jw[n2+jpos] = MIN(jw[n2+jpos], jlev+lev[k]+1);
	    }
	  }
	  else {
	    /* dealing with lower part */
	    if(jpos == -1) {
	      /* this is a fill-in element */
	      jw[lenl] = jcol;
	      jw[n+jcol] = lenl;
	      w[lenl] = -s;
	      jw[n2+lenl] = jlev+lev[k]+1;
	      ++lenl;
	      if(lenl > n) {
		ierr = -1;
		break;
	      }
	    }
	    else {
	      /* this is not a fill-in element */
	      w[jpos] -= s;
	      jw[n2+jpos] = MIN(jw[n2+jpos], jlev+lev[k]+1);
	    }
	  }
	}
      }
      w[jj] = fact;
      jw[jj] = jrow;
    }
    /* reset double-pointer to zero (U-part) */
    if (ii+start < schur_start) {
      for(k = 0; k < lenu; k++) {
	jw[n+jw[ii+start+k]] = -1;
      }
    }
    else {
      for(k = 0; k < lenu; k++) {
	jw[n+jw[schur_start+k]] = -1;
      }

    }

    /* update l-matrix */
    j = 0;
    for(k = 0; k < lenl; k++) {
      if(jw[n2+k] <= lfil) {
	++j;
      }
    }
    data->L->nnzrow[ii+start] = j;
    if(j > 0) {
      // PARMS_NEWARRAY(data->L->pj[ii+start], j);
      // PARMS_NEWARRAY(data->L->pa[ii+start], j);
      data->L->pj[ii+start] = data->L->pj_data + lenl_all;
      data->L->pa[ii+start] = data->L->pa_data + lenl_all;
      lenl_all += j;
      if (lenl_all > max_entr) {
	printf("parms_ilu_vcsr.c, parms_iluk_vcsr: set larger MAX_ROWS_PER_COL and recompile\n");
	MPI_Abort(MPI_COMM_WORLD, 13);
      }
    }
    rowj = data->L->pj[ii+start];
    rowm = data->L->pa[ii+start];
    j = 0;
    for(k = 0; k < lenl; k++) {
      if(jw[n2+k] <= lfil) {
	rowj[j] = jw[k];
	rowm[j] = w[k];
	++j;
      }
    }

    if (ii+start < schur_start) {
/* update u-matrix */
      if(ABS_VALUE(w[ii+start]) <= DBL_EPSILON) {
/*          ierr = ii+1; 
            printf("myid %d: Zero diagonal encountered at row %d \n", rank, ii);
            MPI_Abort(MPI_COMM_WORLD, 13);
*/      
/* Perturb zero diagonal to be (maximum elimination coeff) */
            w[ii+start] = shift;      
      }
      j = 1;
      for(k = ii+1; k < ii+lenu; k++) {
     	  if(jw[n2+start+k] <= lfil) {
	    ++j;
	  }
      }
    }
    else {
      j = 0;
        for(k = schur_start; k < schur_start+lenu; k++) {
	    if(jw[n2+k] <= lfil) {
	      ++j;
	    }
        }
     }

    data->U->nnzrow[ii+start] = j;
    // PARMS_NEWARRAY(data->U->pj[ii+start], j);
    data->U->pj[ii+start] = data->U->pj_data + lenu_all;
    rowj = data->U->pj[ii+start];
    // PARMS_NEWARRAY(data->U->pa[ii+start], j);
    data->U->pa[ii+start] = data->U->pa_data + lenu_all;
    rowm = data->U->pa[ii+start];
    lenu_all += j;

    if (lenu_all > max_entr) {
      printf("parms_ilu_vcsr.c, parms_iluk_vcsr: set larger MAX_ROWS_PER_COL and recompile\n");
      MPI_Abort(MPI_COMM_WORLD, 13);
    }

    PARMS_NEWARRAY(levs[ii+start], j);
    lev = levs[ii+start];
    if (ii+start < schur_start) {
      rowm[0] = 1.0 / w[ii+start];
      rowj[0] = ii+start;
      j = 1;
      for(k = ii+1; k < ii+lenu; k++) {
	if(jw[n2+start+k] <= lfil) {
	  rowj[j] = jw[start+k];
	  rowm[j] = w[start+k];
	  lev[j] = jw[n2+start+k];
	  ++j;
	}
      }
    }
    else {
      j = 0;
      for(k = schur_start; k < schur_start+lenu; k++) {
	if(jw[n2+k] <= lfil) {
	  rowj[j] = jw[k];
	  rowm[j] = w[k];
	  lev[j] = jw[n2+k];
	  ++j;
	}
      }
    }
  }
  for(i = 0; i < m; i++) {
    if (levs[i+start]) {
      PARMS_FREE(levs[i+start]);
    }
  }
  PARMS_FREE(levs);

  data->schur_start = is->schur_start;
  data->n = n;

    if (n-data->schur_start) {
      PARMS_NEWARRAY(data->nnzschur, n-data->schur_start);
      for (i = data->schur_start; i < n; i++) {
	lenl = data->L->nnzrow[i];
	rowj = data->L->pj[i];
	rowm = data->L->pa[i];
	k    = 0;

	for (j = 0; j < lenl; j++) {
	  if (rowj[j] < data->schur_start) {
	    ++k;
	  }
	}
	data->nnzschur[i-data->schur_start] = k;
	ii = 0;
	jj = k;
	for (j = 0; j < lenl; j++) {
	  if (rowj[j] < data->schur_start) {
	    jw[ii]  = rowj[j];
	    w[ii++] = rowm[j];
	  }
	  else {
	    jw[jj]  = rowj[j];
	    w[jj++] = rowm[j];
	  }
	}
	PARMS_MEMCPY(rowj, jw, lenl);
	PARMS_MEMCPY(rowm, w,  lenl);
      }
    }

  PARMS_FREE(w);
  PARMS_FREE(jw);

  /* compute the number of nonzeros in matrix */
  for (i = 0; i < m; i++) {
    data->nnz_mat += mat_vcsr->nnzrow[i];
  }

  /* computer the number of nonzeros in pc */
  for (i = 0; i < m; i++) {
    data->nnz_prec += data->L->nnzrow[i];
    data->nnz_prec += data->U->nnzrow[i];
  }
//  printf("nnzmat = %d, nnzprec = %d \n",data->nnz_mat, data->nnz_prec);
  return ierr;
}

int parms_ilut_vcsr(parms_Mat self, parms_FactParam param, void *mat,
		    parms_Operator *op) 
{
  parms_Operator newOpt;
  parms_ilu_data data;
  parms_vcsr     mat_vcsr;
  parms_Map      is;
  int i, j, j1, j2, ii, jj, n, m, fill, start, schur_start, *jw; 
  int lenl, lenu, len, k, jpos, jrow, jcol, ierr, nnz, *rowj;
  
  REAL dt, tnorm;
  FLOAT *w, *rowm, ft, s, t;

#if defined(DBL_CMPLX)
  double shf = 0.0, ti, sgny;
  int nnzA = 0;    
/*----------get nnz of mat ---------*/
   for(j = 0; j<((parms_vcsr)mat)->n; j++)
	nnzA += ((parms_vcsr)mat)->nnzrow[j];    
#endif  

  is = self->is;
  n = param->n;
  start = param->start;					
  schur_start = param->schur_start;

  if (schur_start == -1) {
    schur_start = n;
  }

  if (!param->isalloc) {
    parms_OperatorCreate(&newOpt);
    PARMS_MEMCPY(newOpt->ops, &parms_ilu_sol_vptr, 1);		
    PARMS_NEW(data);						
    PARMS_NEW(data->L);					
    data->L->n = n;						
    PARMS_NEWARRAY(data->L->nnzrow, n);			
    PARMS_NEWARRAY(data->L->pj, n);				
    PARMS_NEWARRAY(data->L->pa, n);				
    PARMS_NEW(data->U);					
    data->U->n = n;						
    PARMS_NEWARRAY(data->U->nnzrow, n);			
    PARMS_NEWARRAY(data->U->pj,     n);			
    PARMS_NEWARRAY(data->U->pa,     n);
    data->nnz_mat  = 0;
    data->nnz_prec = 0;
    param->isalloc = true;  
    
    newOpt->data = data;
    *op = newOpt;
  }
  else {

    data = (*op)->data;

  }

  mat_vcsr = (parms_vcsr)mat;

  /* malloc memory for working arrays w and jw */
  PARMS_NEWARRAY(w,  2*n);
  PARMS_NEWARRAY(jw, 2*n);

  /* initialize nonzero indicator array */
  for(j = 0; j < n; j++) {
    jw[n+j] = -1;
  }

  fill = param->lfil[0];
  dt   = param->droptol[0];

  /* beginning of main loop */
  ierr = 0;
  m = mat_vcsr->n;
  for (ii = 0; ii < m; ii++) {
    tnorm = 0.0;
    j1 = 0; 
    j2 = mat_vcsr->nnzrow[ii];
    for (k = 0; k < j2; k++) {
      tnorm += ABS_VALUE(mat_vcsr->pa[ii][k]);
    }
    if (ABS_VALUE(tnorm) < 1.0e-20) {
      printf("ilut: zero row encountered at row %d.\n", ii );
      ierr = -3;
      break;
    }
    /* unpack L-part and U-part of row of A in array w */
    lenu = 0;
    lenl = 0;
    if (ii+start < schur_start) {
      jw[ii+start] = ii+start;
      w[ii+start] = 0.0;
      jw[n+ii+start] = ii+start;
      lenu = 1;
      for (j = j1; j < j2; j++) {
	jcol = mat_vcsr->pj[ii][j];
	t = mat_vcsr->pa[ii][j];
	if (jcol < ii+start) {
	  jw[lenl]   = jcol;
	  w[lenl]    = t;
	  jw[n+jcol] = lenl;
	  ++lenl;
	}
	else if(jcol == ii+start){
	  w[ii+start] = t;
	}
	else {
	  jpos       = ii + lenu + start;
	  jw[jpos]   = jcol;
	  w[jpos]    = t;
	  jw[n+jcol] = jpos;
	  ++lenu;
	}
      }
#if defined(DBL_CMPLX)
/*---------------- Add complex shift -- for complex case----------*/
/*------ shift based on droptol -- tau-based shift ----*/	
        shf = dt*tnorm; 
/* ----- shift based on diagonal dominance gap -- dd-based shift -----*/ 
//      shf = schur_start*(tnorm - 2*cabs(w[ii+start]))/(double)nnzA; // rownorm includes diagonal
        ti = cimag(w[ii+start]) ;
        if  (ti<=0) 
	  sgny = -ti - sqrt(ti*ti+shf*shf) ;
       else 
	  sgny = -ti + sqrt(ti*ti+shf*shf) ; 

       w[ii+start] = w[ii+start] + sgny*I;
#endif
       
       tnorm /= (j2-j1);        
    }
    else {
      for (j = j1; j < j2; j++) {
	jcol = mat_vcsr->pj[ii][j];
	t = mat_vcsr->pa[ii][j];
	if (jcol < schur_start) {
	  jw[lenl]   = jcol;
	  w[lenl]    = t;
	  jw[n+jcol] = lenl;
	  ++lenl;
	}
	else {
	  jpos       = schur_start + lenu;
	  jw[jpos]   = jcol;
	  w[jpos]    = t;
	  jw[n+jcol] = jpos;
	  ++lenu;
	}
      }
    }

    len = 0;
    /* eliminate previous rows */
    for (jj = 0; jj < lenl; jj++) {
      /* in order to do the elimination in the correct order we must
	 select the smallest column index among jw[k],
	 k=jj+1,...,lenl.*/
      jrow = jw[jj];
      k    = jj;
      for (j = jj+1; j < lenl; j++) {
	if (jw[j] < jrow) {
	  jrow = jw[j];
	  k = j;
	}
      }
      if (k != jj) {
	/* exchange in jw */
	j = jw[jj];
	jw[jj] = jw[k];
	jw[k] = j;
	/* exchange in jr */
	jw[n+jrow] = jj;
	jw[n+j] = k;
	/* exchange in w */
	s = w[jj];
	w[jj] = w[k];
	w[k] = s;
      }
      jw[n+jrow] = -1;
      /* get the multiplier for row to be eliminated (jrow) */
      rowm = data->U->pa[jrow];
      ft = w[jj]*rowm[0];
      if (ABS_VALUE(ft) <= dt) {
	continue;
      }
      rowj = data->U->pj[jrow];
      nnz  = data->U->nnzrow[jrow];
      /* combine current row and row jrow */
      /* if fill-in element is small then disregard */
      if (ii+start < schur_start) { /* ilut */
	for (k = 1; k < nnz; k++) {
	  s    = ft*rowm[k];
	  jcol = rowj[k];
	  jpos = jw[n+jcol];
	  
	  if (jcol >= ii+start) {
	    /* dealing with upper part */
	    if (jpos == -1) {
	      /* this is a fill-in element */
	      i = ii+start+lenu;
	      jw[i] = jcol;
	      jw[n+jcol] = i;
	      w[i] = -s;
	      ++lenu;
	      if(lenu > n) {
		ierr = -1;
		goto done;
	      }
	    }
	    else {
	      /* this is not a fill-in element */
	      w[jpos] -= s;
	    }
	  }
	  else {
	    /* dealing with lower part */
	    if(jpos == -1) {
	      /* this is a fill-in element */
	      jw[lenl] = jcol;
	      w[lenl] = -s;
	      jw[n+jcol] = lenl;
	      ++lenl;
	      if(lenl > n) {
		ierr = -1;
		goto done;
	      }
	    }
	    else {
	      w[jpos] -= s;
	    }
	  }
	}
      }
      else {		/* partial ilut */
	for (k = 1; k < nnz; k++) {
	  s    = ft*rowm[k];
	  jcol = rowj[k];
	  jpos = jw[n+jcol];
	  
	  if (jcol >= schur_start) {
	    /* dealing with upper part */
	    if (jpos == -1) {
	      /* this is a fill-in element */
	      i = schur_start+lenu;
	      jw[i] = jcol;
	      jw[n+jcol] = i;
	      w[i] = -s;
	      ++lenu;
	      if(lenu > n) {
		ierr = -1;
		goto done;
	      }
	    }
	    else {
	      /* this is not a fill-in element */
	      w[jpos] -= s;
	    }
	  }
	  else {
	    /* dealing with lower part */
	    if(jpos == -1) {
	      /* this is a fill-in element */
	      jw[lenl] = jcol;
	      w[lenl] = -s;
	      jw[n+jcol] = lenl;
	      ++lenl;
	      if(lenl > n) {
		ierr = -1;
		goto done;
	      }
	    }
	    else {
	      w[jpos] -= s;
	    }
	  }
	}
      }

      /* store this pivot element -- (from left to right -- no danger
	  of overlap with the working elements in L (pivots). */
      w[len]  = ft;
      jw[len] = jrow;
      ++len;
    }

    /* reset double-pointer to -1 (U-part) */
    if (ii+start < schur_start) {
      for (k = 0; k < lenu; k++) {
	jw[n+jw[ii+start+k]] = -1;
      }
    }
    else {
      for (k = 0; k < lenu; k++) {
	jw[n+jw[schur_start+k]] = -1;
      }
    }

    /* update l-mat_vcsr */
    lenl = len > fill ? fill : len;
    data->L->nnzrow[ii+start] = lenl;

    /* weigh the elements before sorting */
#if 1
    for (k = 0; k < len; k++) {
      w[k] *= (dt+w[n+jw[k]]);
    }
#endif 
    /* quick sort */
    if (len > lenl) {
      qsplitC(w,jw,len,lenl);
    }
    
    if (lenl > 0) {
      PARMS_NEWARRAY(data->L->pj[start+ii], lenl);
      rowj = data->L->pj[start+ii];
      PARMS_MEMCPY(rowj, jw, lenl);
      PARMS_NEWARRAY(data->L->pa[start+ii], lenl); 
      rowm = data->L->pa[start+ii];
      /* PARMS_MEMCPY(rowm, w, lenl); */
      /*      GCOPY(lenl, w, incx, rowm, incx); */
    }

#if 1
    for (k = 0; k < lenl; k++) {
      rowm[k] = w[k] / (dt + w[n+jw[k]]);
    }
#endif 

    /* update u-part */
    if (ii+start < schur_start) {
      len = 0;
      for (k = 1; k < lenu; k++) {
	 	  if (ABS_VALUE(w[ii+k]) > ABS_VALUE(dt*w[ii])) { 
		    /*
	if (ABS_VALUE(w[ii+start+k]) > dt*tnorm) {
	*/
	  ++len;
	  w[ii+start+len]  = w[ii+start+k];
	  jw[ii+start+len] = jw[ii+start+k];
	}
      }
      lenu = len + 1 > fill ? fill: len + 1;
      jpos = lenu - 1;
      if (len > jpos)
	qsplitC(&w[ii+start+1], &jw[ii+start+1], len, jpos);
      
      data->U->nnzrow[start+ii] = lenu;

      PARMS_NEWARRAY(data->U->pa[start+ii], lenu);
      rowm = data->U->pa[start+ii];

      PARMS_NEWARRAY(data->U->pj[start+ii], lenu);
      rowj = data->U->pj[start+ii];

      /* copy the rest of U */
      if (jpos) {
	PARMS_MEMCPY(&rowj[1], &jw[ii+start+1], jpos);
	PARMS_MEMCPY(&rowm[1], &w[ii+start+1],  jpos);
      }
      t = ABS_VALUE(w[ii+start]);
      for (k = 1; k < lenu; k++) {
	t += ABS_VALUE(w[ii+start+k]);
      }
      w[n+ii+start] = t / (FLOAT)(lenu+1);
      /* store inverse of diagonal element of u */
      if (ABS_VALUE(w[ii+start]) < 1.0e-20) {
	w[ii+start] = tnorm;//(0.0001 + dt)*tnorm;
      }
      rowm[0] = 1.0 / w[ii+start];
      rowj[0] = jw[ii+start];
    }
    else if (ii >= schur_start) {
      data->U->nnzrow[start+ii] = lenu;
      PARMS_NEWARRAY(data->U->pj[start+ii], lenu);
      PARMS_NEWARRAY(data->U->pa[start+ii], lenu);
      PARMS_MEMCPY(data->U->pj[start+ii], &jw[schur_start], lenu);
      PARMS_MEMCPY(data->U->pa[start+ii], &w[schur_start],  lenu);
    }
  }

 done:

  data->schur_start = is->schur_start;
  data->n = n;

    if (n-data->schur_start) {
      PARMS_NEWARRAY(data->nnzschur, n-data->schur_start);
      for (i = data->schur_start; i < n; i++) {
	lenl = data->L->nnzrow[i];
	rowj = data->L->pj[i];
	rowm = data->L->pa[i];
	k    = 0;
	for (j = 0; j < lenl; j++) {
	  if (rowj[j] < data->schur_start) {
	    ++k;
	  }
	}
	data->nnzschur[i-data->schur_start] = k;
	j1 = 0;
	j2 = k;
	for (j = 0; j < lenl; j++) {
	  if (rowj[j] < data->schur_start) {
	    jw[j1]  = rowj[j];
	    w[j1++] = rowm[j];
	  }
	  else {
	    jw[j2]  = rowj[j];
	    w[j2++] = rowm[j];
	  }
	}
	PARMS_MEMCPY(rowj, jw, lenl);
	PARMS_MEMCPY(rowm, w,  lenl);
      }
    }

  PARMS_FREE(w);
  PARMS_FREE(jw);

  /* compute the number of nonzeros in matrix */
  for (i = 0; i < m; i++) {
    data->nnz_mat += mat_vcsr->nnzrow[i];
  }

  /* computer the number of nonzeros in pc */
  for (i = 0; i < m; i++) {
    data->nnz_prec += data->L->nnzrow[i];
    data->nnz_prec += data->U->nnzrow[i];
  }

  return ierr; 
}

/* reuse LU-Facorization (AF) */
int parms_ilu_update(parms_Mat self, parms_FactParam param, void *mat,
		    parms_Operator *op) 
{
  parms_Operator newOpt;
  parms_ilu_data data;
  parms_vcsr     mat_vcsr;
  parms_Map      is;
  int n, start, schur_start, m, iw,*rowjj, nnz, lfil, *jw, **levs, n2, *pj, ierr;
  int lenl, lenu, ii, jj, i, j, k, jcol, jpos, jrow, *rowj, *rowi, *lev, jlev;
  FLOAT *w, *pa, *rowm,*rowmm, t, s, fact, shift, rnrm,t1, t2;
  int rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  is = self->is;
  data = (*op)->data;
  n = data->L->n;
  start = param->start;					
  schur_start = param->schur_start;
  if (schur_start == -1) {
    schur_start = n;
  }
 
  mat_vcsr = (parms_vcsr)mat;

  PARMS_NEWARRAY(jw, n);
  for(i = 0; i < n; i++)
    jw[i] = -1;
  
  for(ii = 0; ii < n; ii++){

    lenl = data->L->nnzrow[ii+start];
    rowj = data->L->pj[ii+start];
    rowm = data->L->pa[ii+start];
    for(j = 0; j < lenl; j++){
      jw[rowj[j]] = j;
      rowm[j] = 0.;
    }
    
    lenu = data->U->nnzrow[ii+start];
    rowjj = data->U->pj[ii+start];
    rowmm = data->U->pa[ii+start];
    if(ii + start < schur_start)
      for(jj = 0; jj < lenu; jj++){
	jw[rowjj[jj]] = jj + ii + start;
	rowmm[jj] = 0.;
      }
    else
      for(jj = 0; jj < lenu; jj++){
	jw[rowjj[jj]] = jj + schur_start;
	rowmm[jj] = 0.;
      }
    
    
    shift = 0.;
    rnrm  = 0.;
    
    nnz = mat_vcsr->nnzrow[ii];    
    pj = mat_vcsr->pj[ii];
    pa = mat_vcsr->pa[ii];

    /*copy row ii in L and U part */
    if(ii + start < schur_start){
      
      for(j = 0; j < nnz; j++){
	rnrm += ABS_VALUE(pa[j]);
	shift = ABS_VALUE(shift) > ABS_VALUE(pa[j]) ? shift : pa[j];
	jcol = jw[pj[j]];    
	if(jcol < 0) 
	  printf("jw %d\n",jcol);	
	else
	  if(jcol < ii + start)
	    rowm[jcol] = pa[j];
	  else
	    rowmm[jcol-ii-start] = pa[j];
      }
      shift /= rnrm;      

    
    }
    else{
      for(j = 0; j < nnz; j++){
	jcol = jw[pj[j]];
	if(jcol < schur_start) 
	  rowm[jcol] = pa[j];
	else
	  rowmm[jcol-schur_start] = pa[j];
      }
    } 	
    

    for(j = 0; j < lenl; j++){
      
      jrow = rowj[j];
      jpos = j;
      
      /* find minimal column index k in row i */
      for(k = j+1; k < lenl; k++)
	if(rowj[k] < jrow){
	  jrow = rowj[k];
	  jpos = k;
	}
      
      if(jpos != j){       
	jcol = rowj[j];
	rowj[j] = rowj[jpos];
	rowj[jpos] = jcol;
	jw[jrow] = j;
	jw[jcol] = jpos;
	t1 = rowm[j];
	rowm[j] = rowm[jpos];
	rowm[jpos] = t1;	
      }      
            
      /* a_{ik} := a_{ik}/a_{kk} */
      nnz   = data->U->nnzrow[jrow];    
      rowjj = data->U->pj[jrow];
      rowmm = data->U->pa[jrow];      
      rowm[j] *= rowmm[0];
      t1 = rowm[j];
      if(ABS_VALUE(t1) < DBL_EPSILON) continue;

      if( ii+start < schur_start ){	
	for(jj = 1; jj < nnz; jj++) {
	  iw = jw[rowjj[jj]];
	  if(iw != -1) 
	    if(iw < ii+start)	      	      
	      rowm[iw] -= t1*rowmm[jj];	  	   
	    else
	      data->U->pa[ii+start][iw-ii-start] -= t1*rowmm[jj];	  
	}
      }
      else {	
	for(jj = 1; jj < nnz; jj++){	  	 
	  iw = jw[rowjj[jj]];
	  if(iw < schur_start)
	    data->L->pa[ii+start][iw] -= t1*rowmm[jj];
	  else
	    data->U->pa[ii+start][iw-schur_start] -= t1*rowmm[jj];
	}	
      }
    }
    if(ii + start < schur_start){
      t1 = data->U->pa[ii+start][0];
      if(ABS_VALUE(t1) < DBL_EPSILON)
	t1 = shift;
      data->U->pa[ii+start][0] = 1./t1;	
    }

    for(i = 0; i < lenl; i++)
      jw[rowj[i]] = -1;
    rowjj = data->U->pj[ii+start];
    for(i = 0; i < lenu; i++)
      jw[rowjj[i]] = -1;     
  }

  
  if(n - data->schur_start){
    PARMS_NEWARRAY0(w, n);
    for(i = data->schur_start; i < n; i++){
      lenl = data->L->nnzrow[i];
      rowj = data->L->pj[i];
      rowm = data->L->pa[i];
      k = data->nnzschur[i-data->schur_start];
      ii = 0; 
      jj = k;
      for(j = 0; j < lenl; j++)
	if(rowj[j] < data->schur_start){
	  jw[ii]  = rowj[j];
	  w[ii++] = rowm[j];
	}
	else{
	  jw[jj]  = rowj[j];
	  w[jj++] = rowm[j];
	}
      PARMS_MEMCPY(rowj, jw,lenl);
      PARMS_MEMCPY(rowm, w, lenl);
    }
    PARMS_FREE(w);    
  }
  
  PARMS_FREE(jw);
  return 0;  
}
