#ifndef _PARMS_SOLVER_IMPL_H_
#define _PARMS_SOLVER_IMPL_H_

#include "parms_operator.h"
#include "parms_solver.h"
#include "parms_mat_impl.h"
#include "parms_pc_impl.h"

typedef struct parms_Solver_ops {
  int (*apply)(parms_Solver self, FLOAT *y, FLOAT *x);
  int (*getresidual)(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *res);
  int (*getresidualnorm2)(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm);
  int (*setksize)(parms_Solver self, int restart);
  int (*setneig)(parms_Solver self, int neig);
  int (*solver_free)(parms_Solver *self);
  int (*solver_view)(parms_Solver self, parms_Viewer v);
} *parms_Solver_ops;

struct parms_Solver_ {

  int ref;
  parms_Solver_ops ops;
  void         *data;
  SOLVERTYPE   stype;
  BOOL         istypeset;
  int          maxits;
  double       tol;
  int          its;
  parms_Mat    A;
  parms_PC     pc;
};

/* external function protos */
extern int fgmres_create(parms_Solver self);
extern int parms_fgmres(parms_Solver self, FLOAT *y, FLOAT *x);

extern int gmres_create(parms_Solver self);
extern int parms_gmres(parms_Solver self, FLOAT *y, FLOAT *x);
/* extern int dgmres_create(parms_Solver self); */

/* additional Solver Type (AF) */
extern int bicgstab_create(parms_Solver self);
extern int parms_bicgstab(parms_Solver self, FLOAT *y, FLOAT *x);

/* additional Solver Type (NR) */
extern int cg_create(parms_Solver self);
extern int parms_cg(parms_Solver self, FLOAT *y, FLOAT *x);
extern int pbicgstab_create(parms_Solver self);
extern int parms_pbicgstab(parms_Solver self, FLOAT *y, FLOAT *x);
extern int pbicgstabras_create(parms_Solver self);
extern int parms_pbicgstabras(parms_Solver self, FLOAT *y, FLOAT *x);
extern int bicgstabras_create(parms_Solver self);
extern int parms_bicgstabras(parms_Solver self, FLOAT *y, FLOAT *x);

#endif 
