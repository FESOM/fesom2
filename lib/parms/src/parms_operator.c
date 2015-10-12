/*--------------------------------------------------------------------

  parms_OperatorCreate 
  parms_OperatorFree 
  parms_OperatorApply 
  parms_OperatorLsol 
  parms_OperatorInvS 
  parms_OperatorAscend 
  parms_OperatorGetSchurPos
  parms_OperatorGetNnz
  parms_OperatorView 

  Users are NOT encouraged to call those functions directly. 

  $Id: parms_operator.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
  ------------------------------------------------------------------*/
#include "include/parms_mat_impl.h"
#include "./include/parms_opt_impl.h"

int parms_OperatorFree(parms_Operator *self)
{

  (*self)->ref--;
  if ((*self)->ref == 0 ) {
    (*self)->ops->operator_free(self);
    PARMS_FREE((*self)->ops);
    PARMS_FREE(*self);
  }
  return 0;
}

int parms_OperatorCreate(parms_Operator *self)
{
  parms_Operator newOpt;

  PARMS_NEW0((newOpt));
  newOpt->ref = 1;
  PARMS_NEW0((newOpt)->ops);

  newOpt->ops->apply  = 0;
  newOpt->ops->lsol   = 0;
  newOpt->ops->invs   = 0;
  newOpt->ops->getu   = 0;
  newOpt->ops->ascend = 0;
  *self = newOpt;
  return 0;
}

int parms_OperatorView(parms_Operator self, parms_Viewer v)
{

  return self->ops->operator_view(self, v);
}

int parms_OperatorApply(parms_Operator self, FLOAT *y, FLOAT *x) 
{
  return self->ops->apply(self, y, x);
}

int parms_OperatorLsol(parms_Operator self, FLOAT *y, FLOAT *x)
{
  return self->ops->lsol(self, y, x);
}
  
int parms_OperatorInvS(parms_Operator self, FLOAT *y, FLOAT *x)
{
  return self->ops->invs(self, y, x);
}

/* used by PCSCHURRAS (AF) */
int parms_OperatorGetU(parms_Operator self, void **mat)
{
  return self->ops->getu(self, mat);
}

int parms_OperatorAscend(parms_Operator self, FLOAT *y, FLOAT *x)
{
  return self->ops->ascend(self, y, x);
}

int parms_OperatorGetSchurPos(parms_Operator self)
{
  return self->ops->getssize(self);
}

void parms_OperatorGetNnz(parms_Operator self, int *nnz_mat, int
			  *nnz_pc)
{
  self->ops->getnnz(self, nnz_mat, nnz_pc);
}
