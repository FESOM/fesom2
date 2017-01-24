#include "include/parms_opt_impl.h"
#include "include/parms_mat_impl.h"

static int MatVec_vcsr(parms_Mat self, FLOAT *x, FLOAT *y)
{
  int          lsize, i, j, *pj, length;
  parms_vcsr   matrix;
  parms_Map    is;
  FLOAT        *pa;
  
  
  /* extract diagonal and off-diagonal matrix */
  matrix   = self->aux_data;
  is       = self->is;
  lsize    = parms_MapGetLocalSize(is);

  for (i = 0; i < lsize; i++) {
    y[i] = 0.0;
    pj  = matrix->pj[i];
    pa  = matrix->pa[i];
    for (j = 0; j <  matrix->nnzrow[i]; j++) {
      y[i]   += pa[j] * x[pj[j]];
    }
  }

  return 0;
}

static int MatMVPY_vcsr(parms_Mat self, FLOAT alpha, FLOAT *x,
			FLOAT beta, FLOAT *y, FLOAT *z) 
{
  int          lsize, i, j, index, *pj, nnz;
  parms_vcsr   matrix;
  parms_Map    is;
  FLOAT        *pa, s;
  
  /* extract diagonal and off-diagonal matrix */
  matrix   = self->aux_data;
  is       = self->is;
  lsize    = is->lsize;

  for (i = 0; i < lsize; i++) {
    s = beta * y[i];
    nnz = matrix->nnzrow[i];
    pj  = matrix->pj[i];
    pa  = matrix->pa[i];
    for (j = 0; j < nnz; j++) {
      index = pj[j];
      s    += alpha * pa[j] * x[index];
    }
    z[i] = s;
  }
  
  return 0;
}

static int MatSetup_vcsr(parms_Mat self)
{
  /* free member space in aux_data */
  if(self->aux_data->space)
    PARMS_FREE(self->aux_data->space);
  self->issetup = true;
  self->is->ninf_send = 0;
  self->is->schur_start =  self->is->lsize;
  return 0;
}

static int MatGetDiag_vcsr(parms_Mat self, void **mat)
{
  if (self->ilutype == PCARMS) {
    parms_vcsr diag, diag_mat;
    int i, nnz;

    diag_mat = self->aux_data;
    PARMS_NEW(diag);
    diag->n = diag_mat->n;
    PARMS_NEWARRAY(diag->nnzrow, diag->n);
    PARMS_NEWARRAY(diag->pj,     diag->n);
    PARMS_NEWARRAY(diag->pa,     diag->n);
    for (i = 0; i < diag->n; i++) {
      nnz = diag->nnzrow[i] = diag_mat->nnzrow[i];
      if (nnz) {
	PARMS_NEWARRAY(diag->pj[i], nnz);
	PARMS_NEWARRAY(diag->pa[i], nnz);
	PARMS_MEMCPY(diag->pj[i], diag_mat->pj[i], nnz);
	PARMS_MEMCPY(diag->pa[i], diag_mat->pa[i], nnz);
      }
    }
    *mat = diag;
    return 0;
  }
  *mat = self->aux_data;
  return 0;
}

static struct parms_Mat_ops parms_matops_vcsr = {
  MatVec_vcsr,
  MatSetup_vcsr,
  0,
  MatMVPY_vcsr,
  MatGetDiag_vcsr,
  MatGetDiag_vcsr,
  0,
  0,
  0,
  0
};

int parms_MatView_vcsr(parms_Mat self, parms_Viewer v)
{
  int        i, j, lsize, pid, length, *rja;
  FLOAT      *ra;
  parms_vcsr matrix;
  parms_Map  is;
  FILE       *fp;

  is      = self->is;
  lsize   = parms_MapGetLocalSize(is);
  pid     = is->pid;
  parms_ViewerGetFP(v, &fp);

  fprintf(fp, "There is one processor available\n");
  fprintf(fp, "The format of output local equations is:\n");
  fprintf(fp, "(pid,local_row_index)=>(pid,global_row_index)\n");
  fprintf(fp, "(pid,local_row_index, local_column_index) = value\n"); 
  fprintf(fp, "\n");

  for (i = 0; i < lsize; i++) {
    fprintf(fp, "(%d,%d) => (%d, %d)\n", pid, i, pid, i);
  }
  
  matrix = self->aux_data;
  fprintf(fp, "Local matrix on processor %d is\n", pid);
#if defined(DBL_CMPLX)  
  for (i = 0; i < lsize; i++) {
    length = matrix->nnzrow[i];
    rja    = matrix->pj[i];
    ra     = matrix->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "(%d,%d,%d) = (%f, %f)  ", pid, i, rja[j], creal(ra[j]), cimag(ra[j]));
    }
    fprintf(fp, "\n");
  }
#else
  for (i = 0; i < lsize; i++) {
    length = matrix->nnzrow[i];
    rja    = matrix->pj[i];
    ra     = matrix->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "(%d,%d,%d) = %f  ", pid, i, rja[j], ra[j]);
    }
    fprintf(fp, "\n");
  }
#endif
  parms_ViewerStoreFP(v, fp);
  return 0;
}

int parms_MatViewCOO_vcsr(parms_Mat self, parms_Viewer v)
{
  int        i, j, lsize, pid, length, *rja;
  FLOAT      *ra;
  parms_vcsr matrix;
  parms_Map  is;
  FILE       *fp;

  is      = self->is;
  lsize   = parms_MapGetLocalSize(is);
  pid     = is->pid;
  parms_ViewerGetFP(v, &fp);

  fprintf(fp, "There is one processor available\n");
  fprintf(fp, "The format of output local equations is:\n");
  fprintf(fp, "(local_row_index local_column_index value\n"); 
  fprintf(fp, "\n");
/*
  for (i = 0; i < lsize; i++) {
    fprintf(fp, "(%d,%d) => (%d, %d)\n", pid, i, pid, i);
  }
*/  
  matrix = self->aux_data;
  fprintf(fp, "Local diagonal matrix on processor %d is\n", pid);
#if defined(DBL_CMPLX)  
  for (i = 0; i < lsize; i++) {
    length = matrix->nnzrow[i];
    rja    = matrix->pj[i];
    ra     = matrix->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "%d  %d  %f  %f  \n", i, rja[j], creal(ra[j]), cimag(ra[j]));
    }
    fprintf(fp, "\n");
  }
#else
  for (i = 0; i < lsize; i++) {
    length = matrix->nnzrow[i];
    rja    = matrix->pj[i];
    ra     = matrix->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "%d  %d  %f  \n", i, rja[j], ra[j]);
    }
    fprintf(fp, "\n");
  }
#endif
  parms_ViewerStoreFP(v, fp);
  return 0;
}


int parms_MatCreate_vcsr(parms_Mat self)
{
  PARMS_MEMCPY(self->ops, &parms_matops_vcsr, 1);
  return 0;
}
