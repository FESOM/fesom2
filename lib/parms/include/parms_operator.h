/**
 * @file   parms_operator.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:51:20 2006
 * 
 * @brief  Functions related to the preconditioner operator objects. 
 * 
 * A operator created implicitly by invoking a ILU-type function. The
 * corresponding destructor is set automatically when the ILU-type
 * function is called.  
 *
 */

#ifndef _PARMS_OPERATOR_H_
#define _PARMS_OPERATOR_H_

#include "parms_sys.h"
#include "parms_vec.h"
#include "parms_viewer.h"

PARMS_CXX_BEGIN

typedef struct parms_Operator_ *parms_Operator;

typedef int (*opt_apply)(parms_Operator self, FLOAT *y, FLOAT *x);

/** 
 * Create an operator object.
 * 
 * @param self A pointer to the operator object created.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorCreate(parms_Operator *self);

/** 
 * Free the memory for the operator object pointed to by self.
 * 
 * @param self A pointer to the operator object.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorFree(parms_Operator *self);

/** 
 * Dump the operator object.
 *
 * Output the L and U part of the operator object.
 * 
 * @param self An operator object.   
 * @param v    A parms_Viewer object.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorView(parms_Operator self, parms_Viewer v);

/** 
 * Perform \f$x = self^{-1}y\f$.
 *
 * Assume \f$ A \approx \left( \begin{array}{ll} L & 0 \\ FU^{-1} &
 * I\end{array} \right) \left( \begin{array}{ll} U & L^{-1}E \\ 0 &
 * S\end{array} \right)\f$, this function performs \f$x =
 * \left(\begin{array}{ll} U & L^{-1}E \\ 0 &
 * S\end{array} \right)^{-1}
 * \left( \begin{array}{ll} L & 0 \\ FU^{-1} & I\end{array}
 * \right)^{-1}y\f$.
 * 
 * @param self An operator.             
 * @param y    A right-hand-side vector.
 * @param x    The solution vector.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorApply(parms_Operator self, FLOAT *y,
			       FLOAT *x); 

/** 
 * Perform forward sweep.
 *
 * Assume \f$ A \approx \left( \begin{array}{ll} L & 0 \\ FU^{-1} &
 * I\end{array} \right) \left( \begin{array}{ll} U & L^{-1}E \\ 0 &
 * S\end{array} \right)\f$, this function performs \f$x =
 * \left( \begin{array}{ll} L & 0 \\ FU^{-1} & I\end{array}
 * \right)^{-1}y\f$.
 * 
 * @param self An operator object.      
 * @param y    A right-hand-side vector.
 * @param x    The solution vector.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorLsol(parms_Operator self, FLOAT *y, FLOAT
			      *x);  

/** 
 * Perform \f$x = S^{-1}y\f$.
 *
 * Assume \f$ A \approx \left( \begin{array}{ll} L & 0 \\ FU^{-1} &
 * I\end{array} \right) \left( \begin{array}{ll} U & L^{-1}E \\ 0 &
 * S\end{array} \right)\f$, this function performs \f$x = S^{-1}y\f$.
 *
 * @param self An operator object.      
 * @param y    A right-hand-side vector.
 * @param x    The solution vector.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorInvS(parms_Operator self, FLOAT *y, FLOAT
			      *x); 



/** 
 * Get local Schur Complement $S_i$ used by PCSCHURRAS (AF)
 *
 * @param self An operator object.      
 * @param mat  local Schur_Complement  
 * 
 * @return 0 on success.
 */
extern int parms_OperatorGetS(parms_Operator self, void **mat); 

/** 
 * Perform block backward substitution.
 *
 * Assume \f$ A \approx \left( \begin{array}{ll} L & 0 \\ FU^{-1} &
 * I\end{array} \right) \left( \begin{array}{ll} U & L^{-1}E \\ 0 &
 * S\end{array} \right)\f$, this function performs \f$x =
 * \left(\begin{array}{ll} U & L^{-1}E \\ 0 &
 * I\end{array} \right)^{-1}y\f$.
 *
 * @param self An operator object.      
 * @param y    A right-hand-side vector.
 * @param x    The solution vector.
 * 
 * @return 0 on success.
 */
extern int parms_OperatorAscend(parms_Operator self, FLOAT *y, FLOAT
				*x); 

/** 
 * Return the start position of the Schur complement in the local
 * matrix. 
 * 
 * @param self An operator object. 
 * 
 * @return The beginning location of the Schur complement in the
 *         local system.
 */
extern int parms_OperatorGetSchurPos(parms_Operator self);

/** 
 * Get the number of nonzero entries of the original matrix and
 * the preconditioning matrix.
 * 
 * @param self     An operator object.                                
 * @param nnz_mat  A pointer to the number of nonzeros of the original
 *                 matrix. 
 * @param nnz_pc   A pointer to the number of nonzeros of the
 *                 preconditioning matrix.    
 */
extern void parms_OperatorGetNnz(parms_Operator self, int *nnz_mat,
				 int *nnz_pc);    

PARMS_CXX_END

#endif
