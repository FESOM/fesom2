#include "parms_solver.h"

#if defined(FORTRAN_CAPS)
#define parms_solverapply_        PARMS_SOLVERAPPLY
#define parms_solvercreate_	  PARMS_SOLVERCREATE
#define parms_solverfree_	  PARMS_SOLVERFREE
#define parms_solvergetits_	  PARMS_SOLVERGETITS
#define parms_solvergetmat_	  PARMS_SOLVERGETMAT
#define parms_solversetparam_	  PARMS_SOLVERSETPARAM
#define parms_solversettype_	  PARMS_SOLVERSETTYPE
#define parms_solvergetpc_      PARMS_SOLVERGETPC
#define parms_solvergetresidual_   PARMS_SOLVERGETRESIDUAL
#define parms_solvergetresidualnorm2_  PARMS_SOLVERGETRESIDUALNORM2
#define parms_solverview_   PARMS_SOLVERVIEW
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_solverapply_        parms_solverapply__
#define parms_solvercreate_	  parms_solvercreate__
#define parms_solverfree_	  parms_solverfree__
#define parms_solvergetits_	  parms_solvergetits__
#define parms_solvergetmat_	  parms_solvergetmat__
#define parms_solversetparam_	  parms_solversetparam__
#define parms_solversettype_	  parms_solversettype__
#define parms_solvergetpc_      parms_solvergetpc__
#define parms_solvergetresidual_   parms_solvergetresidual__
#define parms_solvergetresidualnorm2_  parms_solvergetresidualnorm2__
#define parms_solverview_   parms_solverview__
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_solverapply_        parms_solverapply
#define parms_solvercreate_	  parms_solvercreate
#define parms_solverfree_	  parms_solverfree
#define parms_solvergetits_	  parms_solvergetits
#define parms_solvergetmat_	  parms_solvergetmat
#define parms_solversetparam_	  parms_solversetparam
#define parms_solversettype_	  parms_solversettype
#define parms_solvergetpc_      parms_solvergetpc
#define parms_solvergetresidual_   parms_solvergetresidual
#define parms_solvergetresidualnorm2_  parms_solvergetresidualnorm2
#define parms_solverview_   parms_solverview
#endif

void parms_solverapply_(parms_Solver *self, FLOAT *y, FLOAT
			*x, int *ierr)
{
  *ierr = parms_SolverApply(*self, y, x);
}

void parms_solvercreate_(parms_Solver *self, parms_Mat *A, parms_PC
			 *pc, int *ierr)
{
  *ierr = parms_SolverCreate(self,  *A,  *pc);
}

void parms_solverfree_(parms_Solver *self, int *ierr)
{
  *ierr = parms_SolverFree(self);
}

void parms_solvergetits_(parms_Solver *self, int *its, int *ierr)
{
  *its = parms_SolverGetIts(*self);
  *ierr = 0;
}

void parms_solvergetmat_(parms_Solver *self, parms_Mat *A, int *ierr)
{
  *ierr = parms_SolverGetMat(*self,  A);
}

void parms_solversetparam_(parms_Solver *self, PARAMTYPE *ptype, char
			   *param, int *ierr)
{
  parms_SolverSetParam(*self, *ptype, param);
  *ierr = 0;
}

void parms_solversettype_(parms_Solver *self, SOLVERTYPE *stype, int
			  *ierr)
{
  *ierr = parms_SolverSetType(*self, *stype);
}

void parms_solvergetpc_(parms_Solver *self, parms_PC *PC, int
			  *ierr)
{
  *ierr = parms_SolverGetPC(*self, PC);
}

void parms_solvergetresidual_(parms_Solver *self, FLOAT *y, FLOAT *x, 
                    FLOAT *r, int *ierr)
{
  *ierr = parms_SolverGetResidual(*self, y, x, r);
}

void parms_solvergetresidualnorm2_(parms_Solver *self, FLOAT *y, FLOAT *x, 
                    REAL *rnorm, int *ierr)
{
  *ierr = parms_SolverGetResidualNorm2(*self, y, x, rnorm);
}

void parms_solverview_(parms_Solver *self, parms_Viewer *v, int
			  *ierr)
{
  *ierr = parms_SolverView(*self, *v);
}

