/**
 * @file   parms_pc.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:52:38 2006
 * 
 * @brief  Preconditioner-related Functions.
 * 
 */

#ifndef _PARMS_PC_H_
#define _PARMS_PC_H_

#include "parms_sys.h"
#include "parms_operator.h"
#include "parms_mat.h"
#include "parms_vec.h"
#include "parms_viewer.h"

PARMS_CXX_BEGIN

typedef struct parms_PC_ *parms_PC;
  
/** 
 * Solve \f$self z = y\f$
 * 
 * @param self A preconditioner object. 
 * @param y    A right-hand-side vector.
 * @param z    The solution vector.
 * 
 * @return 0 on success.
 */
extern int parms_PCApply(parms_PC self, FLOAT *y, FLOAT *z);

/** 
 * Set the matrix to create the preconditioning matrix.
 * 
 * @param self A preconditioner object.              
 * @param A    The matrix to be used for creating PC.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetOP(parms_PC self,  parms_Mat A);

/** 
 * Set up the preconditioner (create the preconditioning matrix). 
 * 
 * @param self A preconditioner object.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetup(parms_PC self);

/** 
 * Create a preconditioner object based on the matrix A.
 * 
 * @param self A preconditioner object.
 * @param A    A matrix object.
 * 
 * @return 0 on success.
 */
extern int parms_PCCreate(parms_PC *self, parms_Mat A);

/** 
 * Create an abstract preconditioner object.
 * 
 * @param self A pointer to the preconditioner object.
 * 
 * @return 0 on success.
 */
extern int parms_PCCreateAbstract(parms_PC *self);

/** 
 * Free the memory for the preconditioner object pointed to by self.
 * 
 * @param self A pointer to the memory for the preconditioner object.
 * 
 * @return 0 on success.
 */
extern int parms_PCFree(parms_PC *self);

/** 
 * Dump preconditioner object self.
 * 
 * @param self A preconditioner object.
 * @param v    A viewer object.
 */
extern void parms_PCView(parms_PC self, parms_Viewer v);

/** 
 * Set preconditioner type.
 *
 *  Currently supported preconditioners:
 *  \f{tabular}{ll}
 *  PCBJ     & block Jacobi \\
 *  PCRAS    & restricted additive Schwarz \\
 *  PCSCHUR  & Schur complement 
 *  \f}
 * 
 * @param self   A preconditioner object.         
 * @param pctype The type of preconditioner.      
 *               - PCBJ      block Jacobi 
 *               - PCRAS     restricted additive Schwarz 
 *               - PCSCHUR   Schur complement 
 * 
 * @return 0 on success.
 */
extern int parms_PCSetType(parms_PC self, PCTYPE pctype);

/** 
 * Set local preconditioner type.
 *
 * Supported ILU preconditioners:
 * 
 *  \f{tabular}{ll}
 *  PCILU0    &  ILU0  \\
 *  PCILUK    &  ILUK  \\
 *  PCILUT    &  ILUT in SPARSKIT \\
 *  PCARMS    &  ARMS implemented by Yousef Saad
 *  \f}
 
 * @param self     A preconditioner object.         
 * @param pcstype  The type of local preconditioner:
 *                  - PCILU0
 *                  - PCILUK
 *                  - PCILUT
 *                  - PCARMS
 *
 * @return 0 on success.
 */
extern int parms_PCSetILUType(parms_PC self, PCILUTYPE pcstype);

/** 
 * Perform ILU factorization.
 * 
 * @param self   A preconditioner object.                                  
 * @param param  The parameters used for ILU factorization.        
 * @param mat    Matrix to be factored. 
 * @param op     The operator created.   
 * 
 * @return 0 on success.
 */
extern int parms_PCILU(parms_PC self, parms_FactParam param, void
			*mat, parms_Operator *op);

/** 
 * Set parameters for the preconditioner object.
 *
 * Supported parameters:
 *   - tol    drop tolerance
 *   - fil    fill-in
 *   - nlev   number of levels
 *   - bsize  block size for finding independent sets in ARMS.
 *   - tolind drop tolerance for finding independent sets.
 *   - iksize the restart size for the inner GMRES.
 *   - imax   the number of iterations for the inner GMRES.
 *
 * @param self   A preconditioner object.
 * @param nflags The number of parameters.
 * @param params A pointer to parameters.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetParams(parms_PC self, int nflags, char **params);

/** 
 * Set permutation and scaling options for interlevel blocks.
 * 
 * @param self A preconditioner object.
 * @param meth Options:
 *             - meth[0] nonsummetric permutations of  1: yes. affects
 *                       rperm USED FOR LAST SCHUR COMPLEMENT 
 *             - meth[1] permutations of columns 0:no 1: yes. So far
 *                       this is USED ONLY FOR LAST BLOCK [ILUTP
 *                       instead of ILUT]. (so ipar[11] does not matter
 *                       ,enter zero). If ipar[15] is one then ILUTP
 *                       will be used instead of ILUT. Permutation data
 *                       stored in: perm2.  
 *              - meth[2] diag. row scaling. 0:no 1:yes. Data: D1
 *              - meth[3] diag. column scaling. 0:no 1:yes. Data: D2
 * @param flag Options:
 *              - 1 interlevel block.
 *              - 0 last block.
 *
 * @return 0 on success.
 */
extern int parms_PCSetPermScalOptions(parms_PC self, int *meth, int
				      flag); 

/** 
 * Set fill-in parameter for ILUT and ARMS.
 * 
 * @param self A preconditioner object.                
 * @param fill A int array of size 7.
 *             - fill[0] amount of fill-in kept  in \f$L_{B}\f$. 
 *             - fill[1] amount of fill-in kept  in \f$U_{B}\f$.
 *             - fill[2] amount of fill-in kept  in \f$E L^{-1}\f$.
 *             - fill[3] amount of fill-in kept  in \f$U^{-1}_{B} F\f$.
 *             - fill[4] amount of fill-in kept  in \f$S\f$.
 *             - fill[5] amount of fill-in kept  in \f$L_S\f$.
 *             - fill[6] amount of fill-in kept  in \f$U_S\f$. 
 * 
 * @return 0 on success.
 */
extern int parms_PCSetFill(parms_PC self, int *fill);
  
/** 
 * Set the drop tolerance for ILUT preconditioner.
 * 
 * @param self A preconditioner object.
 * @param tol  A double array of size 7. 
 *             - tol[0]  threshold for dropping  in L_{B}. 
 *             - tol[1]  threshold for dropping  in U_{B}.
 *             - tol[2]  threshold for dropping  in L^{-1} F 
 *             - tol[3]  threshold for dropping  in E U^{-1} 
 *             - tol[4]  threshold for dropping  in Schur complement
 *             - tol[5]  threshold for dropping  in L in last block
 *             - tol[6]  threshold for dropping  in U in last block
 * 
 * @return 0 on success.
 */
extern int parms_PCSetTol(parms_PC self, double *tol);

/** 
 * Set the number of levels for ILUK and ARMS.
 * 
 * @param self   A preconditioner object.
 * @param nlevel The number of levels.     
 * 
 * @return 0 on success.
 */
extern int parms_PCSetNlevels(parms_PC self, int nlevel);

/** 
 * Set the type of permutation in ARMS
 * 
 * @param self A preconditioner object.
 * @param type Permutation type.
 *             - 1 non-symmetric permutaion.
 *             - 0 symmetric permutation.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetPermType(parms_PC self, int type);

/** 
 * Set the block size for ARMS.
 * 
 * @param self  A preconditioner object.
 * @param bsize The block size for ARMS.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetBsize(parms_PC self, int bsize);

/** 
 * Set the restart size for the inner GMRES.
 * 
 * @param self A preconditioner object.
 * @param im   The restart size of the inner GMRES.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetInnerKSize(parms_PC self, int im);


/** 
 * Set the maximum iteration counts for the inner GMRES.
 * 
 * @param self A preconditioner object.     
 * @param imax The maximum iteration counts.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetInnerMaxits(parms_PC self, int imax);

/** 
 * Set the convergence tolerance for the inner GMRES.
 * 
 * @param self A preconditioner object.  
 * @param eps  The convergence tolerance.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetInnerEps(parms_PC self, REAL eps);

/** 
 * Set the tolerance for finding independent sets.
 * 
 * @param self    A preconditioner object.
 * @param tolind  The drop tolerance for finding independent sets.
 * 
 * @return 0 on success.
 */
extern int parms_PCSetTolInd(parms_PC self, REAL tolind);

/** 
 * Get the ratio of the number of nonzero entries of the
 * preconditioning matrix to that of the original matrix.
 * 
 * @param self  A preconditioner.            
 * @param ratio A pointer to the ratio.      
 * 
 * @return 0 on success.
 */
extern int parms_PCGetRatio(parms_PC self, double *ratio);

/** 
 * Return the name of a preconditioner.
 * 
 * @param self A preconditioner.               
 * @param name The name of preconditioner.     
 * 
 * @return 0 on success.
 */
extern int parms_PCGetName(parms_PC self, char **name);

/** 
 * Return the name of a local preconditioner.
 * 
 * @param self    A preconditioner.                         
 * @param iluname The name of local ILU preconditioner.     
 * 
 * @return 0 on success.
 */
extern int parms_PCILUGetName(parms_PC self, char **iluname);

/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_pccreate_(parms_PC *self, parms_Mat *A, int *ierr);

extern void parms_pcfree_(parms_PC *self, int *ierr);

extern void parms_pcgetratio_(parms_PC *self, double *ratio, int *ierr);

extern void parms_pcsetbsize_(parms_PC *self, int *bsize, int *ierr);
               
extern void parms_pcsetfill_(parms_PC *self, int *fill, int *ierr);

extern void parms_pcsetilutype_(parms_PC *self, PCILUTYPE *pcstype, int
			 *ierr);
			 
extern void parms_pcsetinnereps_(parms_PC *self, REAL *eps, int *ierr);

extern void parms_pcsetinnerksize_(parms_PC *self, int *im, int *ierr);

extern void parms_pcsetinnermaxits_(parms_PC *self, int *imax, int *ierr);

extern void parms_pcsetnlevels_(parms_PC *self, int *nlevel, int *ierr);

extern void parms_pcsetparams_(parms_PC *self, int *nflags, char **params,
			int *ierr);

extern void parms_pcsettol_(parms_PC *self, double *tol, int *ierr);

extern void parms_pcsettolind_(parms_PC *self, REAL *tolind, int *ierr);

extern void parms_pcsettype_(parms_PC *self, PCTYPE *pctype, int *ierr);

extern void parms_pcsetup_(parms_PC *self, int *ierr);

extern void parms_pcsolve_(parms_PC *self, FLOAT *y, FLOAT *z, int *ierr);

extern void parms_pcview_(parms_PC *self, parms_Viewer *v, int *ierr);

extern void parms_pcgetname_(parms_PC *self, char *name, int *size, int
		      *ierr, int len);
		      
extern void parms_pcilugetname_(parms_PC *self, char *name, int *size, int
			 *ierr, int len);
			 
extern void parms_pccreateabstract_(parms_PC *self, int *ierr);

extern void parms_pcsetpermscaloptions_(parms_PC *self, int *meth, int
				      *flag, int *ierr);

extern void parms_pcsetpermtype_(parms_PC *self, int *type, int *ierr);

extern void parms_pcsetop_(parms_PC *self, parms_Mat *A, int *ierr);

/*
 *
 * end Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif 
