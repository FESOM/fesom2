#include "./include/parms_pc_impl.h"
#include "./include/parms_opt_impl.h"

typedef struct bj_data {
  parms_Operator op;
  BOOL issetup;
} *bj_data;

/** Free the memory for struct bj_data.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
static int pc_bj_free(parms_PC *self)
{
  bj_data pc_data;
  parms_Operator op;

  pc_data = (bj_data)(*self)->data;
  op = pc_data->op;
  parms_OperatorFree(&op);
  PARMS_FREE(pc_data);
  (*self)->param->isalloc = false;  
  return 0;
}

/** Dump BJ preconditioner.
 *
 *  \param self  A preconditioner object.
 *  \param v     A viewer object.
 *  \return 0 on success.
 */
static int pc_bj_view(parms_PC self, parms_Viewer v)
{
  bj_data pc_data;

  pc_data = (bj_data)self->data;
  parms_OperatorView(pc_data->op, v);
  return 0;
}

/** Set up block Jacobi preconditioner.
 *
 *  \param self  A preconditioner object.
 *  \return 0 on success.
 */
static int pc_bj_setup(parms_PC self)
{
  bj_data pc_data;
  parms_Mat A;
  void *diag_mat;
  parms_FactParam param;
  parms_Operator op;

  A =  self->A;
  param = self->param; 
  pc_data = (bj_data)self->data;
  if(pc_data->issetup)
    op = pc_data->op;
  else{
    param->start = 0;
    param->n = A->is->lsize;
    param->schur_start = -1;
    pc_data->issetup = true;
    op = NULL;
  }

  /* reuse lu-factorization (AF) */
  parms_MatGetDiag(A, &diag_mat);
  parms_PCILU(self, param, diag_mat, &op);


  pc_data->op = op;
 
  return 0;
}

/** Apply preconditioner BJ to the vector y.
 *
 *  \x = self^{-1}y.
 *
 *  \param self  A preconditioner object.
 *  \param y     A right-hand-side vector.
 *  \param x     Solution vector.
 *  \return 0
 */
static int pc_bj_apply(parms_PC self, FLOAT *y, FLOAT *x)
{
  bj_data pc_data;
  parms_Operator op;

  pc_data = (bj_data)self->data;
  op = pc_data->op;
  parms_OperatorApply(op, y, x);
  return 0;
}

/** Get the ratio of the number of nonzero entries of the
 *  preconditioning matrix to that of the original matrix.
 *
 *  \param self   A preconditioner.
 *  \param ratio  A pointer to the ratio.      
 *  \return 0 on success.
 */
static int pc_bj_getratio(parms_PC self, double *ratio)
{
  bj_data pc_data;
  parms_Operator op;
  int nnz_mat, nnz_pc;
  int gnnz_mat, gnnz_pc;

  pc_data = (bj_data)self->data;
  op = pc_data->op;
  parms_OperatorGetNnz(op, &nnz_mat, &nnz_pc);
  MPI_Allreduce(&nnz_mat, &gnnz_mat, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD); 
  MPI_Allreduce(&nnz_pc, &gnnz_pc, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD); 
  *ratio = (double)gnnz_pc/(double)gnnz_mat;  
  return 0;
}

/** Create a block Jacobi preconditioner.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success
 */
int parms_PCCreate_BJ(parms_PC self)
{
  bj_data data;

  PARMS_NEW(data);
  data->issetup = false;

  self->data = data;
  self->ops->pc_free     = pc_bj_free;
  self->ops->pc_view     = pc_bj_view;
  self->ops->apply    = pc_bj_apply;
  self->ops->setup    = pc_bj_setup;
  self->ops->getratio = pc_bj_getratio;
  return 0;
}

