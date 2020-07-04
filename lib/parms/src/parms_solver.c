/*--------------------------------------------------------------------
  parms_SolverApply    : solve the linear system of equation 
  parms_SolverCreate   : create a parms_Solver object.
  parms_SolverFree     : free the memory for the parms_Solver object.
  parms_SolverGetIts   : get the iteration counts.
  parms_SolverGetMat   : get the matrix for the linear system solved.
  parms_SolverGetPC    : get the preconditioning matrix.
  parms_SolverSetParam : set the parameters for the solver (maxits,
                         tolerance, etc.)
  parms_SolverSetType  : set the type of the solver. FGMRES and DGMRES
                         are supplied in the package.
                         dgmres removed YS 
  parms_SolverView     : dump the maximum iteration count and the
                         iteration counts

  A code fragment for using solver functions:

  parms_Solver solver;
  parms_Mat A;
  parms_PC  pc;

  // create a solver 
  parms_SolverCreate(&solver, A, pc);

  // set the maximum number of iterations  
  parms_SolverSetParam(solver, MAXITS, "300");
  // set the restart size of FGMRES or DGMRES
  parms_SolverSetParam(solver, KSIZE, "60");
  // set the convergence tolerance 
  parms_SolverSetParam(solver, DTOL, "1.0e-6");
  // set the number of eigenvectors
  parms_SolverSetParam(solver, NEIG, "8"); 

  // solver the linear system 
  parms_SolverApply(solver, rhs, x);

  // free the memory for solver
  parms_SolverFree(&solver);

  $Id: parms_solver.c,v 1.5 2006-12-01 20:44:20 zzli Exp $
 -------------------------------------------------------------------*/
#include <stdlib.h>
#include "parms_vec.h"
#include "parms_viewer.h"
#include "parms_mat_impl.h"
#include "parms_pc_impl.h"
#include "parms_solver_impl.h"

/*int fgmres_create(parms_Solver self);*/
/* extern int dgmres_create(parms_Solver self); */

int parms_SolverView(parms_Solver self, parms_Viewer v)
{
  self->ops->solver_view(self, v);

  return 0;
}

/** 
 * Create a parms_Solver object.
 * 
 * @param self A pointer to the parms_Solver object created.
 * @param A    The matrix of the linear system.   
 * @param pc   The preconditioner.  
 * 
 * @return 0 on success.
 */

int parms_SolverCreate(parms_Solver *self, parms_Mat A, parms_PC pc)
{
  parms_Solver new_solver;

  PARMS_NEW0((new_solver));
  new_solver->ref = 1;
  PARMS_NEW0((new_solver)->ops);
  new_solver->ops->apply    = 0;
  new_solver->ops->setksize = 0;
  new_solver->istypeset = false;
  new_solver->A = A;
  A->ref++;
  new_solver->pc = pc;
  pc->ref++;
  new_solver->maxits = 100;
  new_solver->tol    = 1.0e-6;
  *self = new_solver;
  return 0;
}

/** 
 * Solve the equation \f$Ax = y\f$.
 * 
 * @param self A parms_Solver object.     
 * @param x    The solution vector.       
 * @param y    The right-hand-side vector.
 * 
 * @return 0 on success.
 */
int parms_SolverApply(parms_Solver self, FLOAT *y, FLOAT *x)
{
//  if(self->istypeset == true){
//      self->ops->solver_free(&self);
//      self->istypeset = false;
//  }
  if (self->istypeset == false) {
/* Default solver - fgmres */
	parms_SolverSetType(self, SOLFGMRES);
  }

  return self->ops->apply(self, y, x);
}

/** 
 * Compute the local residual \f$ r = y - Ax \f$.
 * 
 * @param self A parms_Solver object.     
 * @param x    The solution vector.       
 * @param y    The right-hand-side vector.
 * @param r    The computed residual vector.
 * 
 * @return 0 on success.
 */

int parms_SolverGetResidual(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *r)
{
  return self->ops->getresidual(self, y, x, r);
}

/** 
 * Compute the local residual 2-norm \f$ ||r = y - Ax|| \f$.
 * 
 * @param self A parms_Solver object.     
 * @param x    The solution vector.       
 * @param y    The right-hand-side vector.
 * @param rnorm The 2-norm of the residual vector.
 * 
 * @return 0 on success.
 */

int parms_SolverGetResidualNorm2(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm)
{
  return self->ops->getresidualnorm2(self, y, x, rnorm);
}
  
/** 
 * Set the type of the solver.
 *
 * Only FGMRES solver is available in the package.
 * 
 * @param self   A parms_Solver object.       
 * @param stype  The type of Krylov subspace. 
 *               -SOLFGMRES
 *               -SOLDGMRES
 *               
 * @return 0 on success.
 */

int parms_SolverSetType(parms_Solver self, SOLVERTYPE stype)
{
  if (self->istypeset && self->stype == stype) {
    return 0;
  }
  if (self->istypeset) {
    self->ops->solver_free(&self);
  }
  if(stype == SOLFGMRES)
	fgmres_create(self);
  else if(stype == SOLGMRES)
	gmres_create(self);
  else if(stype == SOLBICGS)
        bicgstab_create(self);
  else if(stype == SOLPBICGS)
        pbicgstab_create(self);
  else if(stype == SOLPBICGS_RAS)
        pbicgstabras_create(self);
  else if(stype == SOLBICGS_RAS)
        bicgstabras_create(self);
  else if(stype == SOLCG)
        cg_create(self);
  else{
	printf("ERROR: Invalid choice of solver - (Check SOLVERTYPE for parms_SolverSetType(...) \n");
	PARMS_ABORT(17);
   }
  self->stype      = stype;
  self->istypeset  = true;
  return 0;
}
 
/** 
 * Set parameter for the solver.
 * 
 * Set the maximum iteration counts, the restart size of GMRES, and
 * the convergence tolerance.
 * 
 * @param self   A parms_Solver object.
 * @param ptype  The type of parameter.
 *               -MAXITS maximum iteration counts.
 *               -KSIZE  restart size of GMRES.
 *               -DTOL   converence tolerance.
 *               -NEIG   number of eigenvectors.
 * @param param  Parameters for the solver.
 */
void parms_SolverSetParam(parms_Solver self, PARAMTYPE paramtype, char
			  *param)  
{

  if (self->istypeset == false) {
    parms_SolverSetType(self, SOLFGMRES);
  }
  if (paramtype == MAXITS) {
    self->maxits = atoi(param);
  }
  else if (paramtype == KSIZE) {
    self->ops->setksize(self, atoi(param));
  }
  else if (paramtype == DTOL) {
    self->tol = strtod(param, (char **)NULL);
  }
  else if (paramtype == NEIG) {
    self->ops->setneig(self, atoi(param));
  }
}

/** 
 * Free the memory for the parms_Solver object.
 * 
 * @param self A pointer to the parms_Solver object to be freed. 
 * 
 * @return 0 on success.
 */
int parms_SolverFree(parms_Solver *self)
{
  (*self)->ref--;
  if ((*self)->ref == 0 ) {
    parms_MatFree(&(*self)->A);
    parms_PCFree(&(*self)->pc);
    (*self)->ops->solver_free(self);
    PARMS_FREE((*self)->ops);
    PARMS_FREE(*self);
  }
  return 0;
}

/** 
 * Get the iteration counts.
 * 
 * @param self A parms_Solver object.
 * 
 * @return The iteration counts.
 */
int parms_SolverGetIts(parms_Solver self)
{
  return self->its;
}

/** 
 * Get the matrix of the linear system.
 * 
 * @param self A parms_Solver object.              
 * @param A    A pointer to the matrix returned.   
 * 
 * @return 0 on success.
 */
int parms_SolverGetMat(parms_Solver self, parms_Mat *A)
{

  if (self->ref == 1 ) {
    *A = self->A;
    return 0;
  }
  else {
    *A = NULL;
    return -1;
  }
}

/** 
 * Get the preconditioning matrix.
 * 
 * @param self A parms_Solver object.                    
 * @param PC   A pointer to the preconditioning matrix.  
 * 
 * @return 0 on success.
 */
int parms_SolverGetPC(parms_Solver self, parms_PC *PC)
{
  if (self->ref == 1 ) {
    *PC = self->pc;
    return 0;
  }
  else {
    *PC = NULL;
    return -1;
  }
}

