/**
 * @file   parms_solver.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:58:22 2006
 * 
 * @brief  Functions related to the Krylov subspace methods. Only
 *         FGMRES is supported. 
 * 
 */

#ifndef _PARMS_SOLVER_H_
#define _PARMS_SOLVER_H_

#include "parms_vec.h"
#include "parms_mat.h"
#include "parms_viewer.h"
#include "parms_pc.h"
#include "parms_operator.h"

PARMS_CXX_BEGIN

typedef struct parms_Solver_ *parms_Solver;

/** 
 * Get the matrix of the linear system.
 * 
 * @param self A parms_Solver object.              
 * @param A    A pointer to the matrix returned.   
 * 
 * @return 0 on success.
 */
extern int parms_SolverGetMat(parms_Solver self, parms_Mat *A);

/** 
 * Get the preconditioning matrix.
 * 
 * @param self A parms_Solver object.                    
 * @param PC   A pointer to the preconditioning matrix.  
 * 
 * @return 0 on success.
 */
extern int parms_SolverGetPC(parms_Solver self, parms_PC *PC);

/** 
 * Solve the equation \f$Ax = y\f$.
 * 
 * @param self A parms_Solver object.     
 * @param x    The solution vector.       
 * @param y    The right-hand-side vector.
 * 
 * @return 0 on success.
 */
extern int parms_SolverApply(parms_Solver self, FLOAT *x, FLOAT *y); 

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
extern int parms_SolverGetResidual(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *r);

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
extern int parms_SolverGetResidualNorm2(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm);


/** 
 * Set the type of the solver.
 *
 * Only FGMRES solver is available in the package.
 * 
 * @param self   A parms_Solver object.       
 * @param stype  The type of Krylov subspace. 
 *               - SOLFGMRES
 *               - SOLDGMRES
 *               
 * @return 0 on success.
 */
extern int parms_SolverSetType(parms_Solver self, SOLVERTYPE stype);

/** 
 * Create a parms_Solver object.
 * 
 * @param self A pointer to the parms_Solver object created.
 * @param A    The matrix of the linear system.   
 * @param pc   The preconditioner.  
 * 
 * @return 0 on success.
 */
extern int parms_SolverCreate(parms_Solver *self, parms_Mat A, parms_PC pc);

/** 
 * Free the memory for the parms_Solver object.
 * 
 * @param self A pointer to the parms_Solver object to be freed. 
 * 
 * @return 0 on success.
 */
extern int parms_SolverFree(parms_Solver *self);

/** 
 * Dump te solver object.
 * 
 * @param self A parms_Solver object.   
 * @param v    A parms_Viewer object.   
 */
extern int parms_SolverView(parms_Solver self, parms_Viewer v);

/** 
 * Set parameter for the solver.
 * 
 * Set the maximum iteration counts, the restart size of GMRES, and
 * the convergence tolerance.
 * 
 * @param self   A parms_Solver object.
 * @param ptype  The type of parameter.
 *               - MAXITS maximum iteration counts.
 *               - KSIZE  restart size of GMRES.
 *               - DTOL   converence tolerance.
 *               - NEIG   number of eigenvectors.
 * @param param  Parameters for the solver.
 */
extern void parms_SolverSetParam(parms_Solver self, PARAMTYPE ptype,
				 char *param);
   
/** 
 * Get the iteration counts.
 * 
 * @param self A parms_Solver object.
 * 
 * @return The iteration counts.
 */
extern int parms_SolverGetIts(parms_Solver self);


/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_solverapply_(parms_Solver *self, FLOAT *x, FLOAT
			*y, int *ierr);

extern void parms_solvercreate_(parms_Solver *self, parms_Mat *A, parms_PC
			 *pc, int *ierr);

extern void parms_solverfree_(parms_Solver *self, int *ierr);

extern void parms_solvergetits_(parms_Solver *self, int *its, int *ierr);

extern void parms_solvergetmat_(parms_Solver *self, parms_Mat *A, int *ierr);

extern void parms_solversetparam_(parms_Solver *self, PARAMTYPE *ptype, char
			   *param, int *ierr);

extern void parms_solversettype_(parms_Solver *self, SOLVERTYPE *stype, int
			  *ierr);
			  
extern void parms_solvergetpc_(parms_Solver *self, parms_PC *PC, int
			  *ierr);

extern void parms_solvergetresidual_(parms_Solver *self, FLOAT *y, FLOAT *x, 
                    FLOAT *r, int *ierr);

extern void parms_solvergetresidualnorm2_(parms_Solver *self, FLOAT *y, FLOAT *x, 
                    REAL *rnorm, int *ierr);

extern void parms_solverview_(parms_Solver *self, parms_Viewer *v, int
			  *ierr);

/*
 *
 * End Fortran Wrapper Functions 
 *
*/


PARMS_CXX_END

#endif 
