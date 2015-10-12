/**
 * @file   parms_mat.h
 * @author Zhongze Li
 * @date   Tue Oct 17 10:05:30 2006
 * 
 * @brief  Functions related to the matrix computations.
 * 
 * 
 */

#ifndef _PARMS_MAT_H_
#define _PARMS_MAT_H_

#include "parms_sys.h"
#include "parms_vec.h"
#include "parms_viewer.h"
#include "parms_operator.h"

PARMS_CXX_BEGIN

typedef struct parms_Mat_ *parms_Mat;

/** 
 * Create a parms_Mat object.
 *
 * Create a parms_Mat object based on data distribution layout map.
 * 
 * @param self  A pointer to the parms_Mat object created. 
 * @param map   A parms_Map object, which describes the data
 *              distribution among processors.  
 * 
 * @return 0 on success.
 */
extern int parms_MatCreate(parms_Mat *self, parms_Map map);

/** 
 * Dump parms_Mat object.
 * 
 * @param self A pointer to a parms_Mat object.
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
extern int parms_MatView(parms_Mat self, parms_Viewer v);

extern int parms_MatViewCOO(parms_Mat self, parms_Viewer v);

/** 
 * Free the parms_Mat object pointed to by self.
 * 
 * @param self A pointer to a parms_Mat object.
 * 
 * @return 0 on success.
 */
extern int parms_MatFree(parms_Mat *self);

/** 
 * Perform \f$y = self \times x\f$.
 * 
 * @param self A parms_Mat object.      
 * @param x    A Vector object.      
 * @param y    Another Vector object.
 * 
 * @return 0 on success.
 */
extern int parms_MatVec(parms_Mat self, FLOAT *x, FLOAT *y);

/** 
 * Set up parms_Mat object self.
 *
 * This is the most important functoin for the parms_Mat object. This
 * function combines the function bdry and setup in the old version 
 * of pARMS. The function sets up the data structure needed by the
 * distributed matrix-vector multiplication, divides the variables on
 * the local processors into two categories: interior and interface
 * variables.
 * 
 * @param self A parms_Mat object. 
 * 
 * @return 0 on success.
 */
extern int parms_MatSetup(parms_Mat self);

/** 
 * Insert/add values to the parms_Mat object self.
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	  The number of rows inserted.    
 * @param im 	An array of row global indices.                
 * @param ia 	  An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                    
 * @param mode 	  Insert value mode:                     
 * 		  - INSERT  insert values to parm_Mat self.
 * 		  - ADD     add values to parm_Mat self.
 * 		  
 * @return 0 on success.
 */
extern int parms_MatSetValues(parms_Mat self, int m, int *im, int *ia,
			      int *ja, FLOAT *values, INSERTMODE mode);

/** 
 * Insert/add values to the parms_Mat object self. This assumes matrix 
 * values are being added element-by-element
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	The number of rows inserted.      
 * @param im 	An array of row global indices.                                                 
 * @param ia 	An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                    
 * @param mode 	Insert value mode:                     
 * 		      - INSERT  insert values to parm_Mat self.
 * 		      - ADD     add values to parm_Mat self.
 *
 * NOTE: New entries are always inserted first, so mode does not
 *       matter if this is a new entry. Subsequent calls will either 
 *       replace (mode = INSERT) or add (mode = ADD) to the existing 
 *       entries in a particular position 
 * 		  
 * @return 0 on success.
 */
extern int parms_MatSetElementMatrix(parms_Mat self, int m, int *im, int *ia,
		       int *ja, FLOAT *values, INSERTMODE mode); 
		       
/** 
 * Assembles the finite element matrix by updating off-processor 
 * contributions.
 * 
 * @param self    A parms_Mat object.                                           
 * @return 0 on success.
 */
extern int parms_MatAssembleElementMatrix(parms_Mat self);		       

/** 
 * Insert values to the parms_Mat object self. This assumes matrix 
 * values for this row have already been set, and are to be replaced 
 * by the new ones provided as input.
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	The number of rows inserted.      
 * @param im 	An array of row global indices.                                                 
 * @param ia 	An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                      
 * @return 0 on success.
 */
extern int parms_MatResetRowValues(parms_Mat self, int m, int *im, int *ia,
		       int *ja, FLOAT *values);

/** 
 * Reset the matrix to be re-used. 
 * @param self    A parms_Mat object.                                
 * @param nonzerostructure  The nonzero structure:
 *                        SAME_NONZERO_STRUCTURE
 *                        DIFFERENT_NONZERO_STRUCTURE     
 *                     
 * @return 0 on success.
 */
extern int parms_MatReset(parms_Mat self, NNZSTRUCT nonzerostructure);

/** 
 * Set the communication type.
 *
 * Set the communication style across processors.
 * communication style:
 *  - P2P       point-to-point (data copied to/from auxilliary buffers).
 *  - DERIVED   derived datatype.
 *
 * @param self  A matrix object.
 * @param ctype Communication style:
 *              - P2P     point-to-point (data copied to/from
 *                        auxilliary buffers).
 *              - DERIVED derived datatype.
 *              
 * @return 0 on success.
 */
extern int parms_MatSetCommType(parms_Mat self, COMMTYPE ctype);

/** 
 * Get the diagonal part of the local matrix. 
 * 
 * @param self A parms_Mat object.
 * @param mat  The diagonal part of the local matrix.    
 * 
 * @return 0 on success.
 */
extern int parms_MatGetDiag(parms_Mat self, void **mat);


/** 
 * Perform \f$z = alpha*self*x + beta*y\f$.
 * 
 * @param self   A matrix object.           
 * @param alpha  A scalar.                  
 * @param x 	 A vector object.           
 * @param beta 	 A scalar.                  
 * @param y 	 A vector object.           
 * @param z 	 A vector stores the result.
 * 
 * @return 0 on success.
 */
extern int parms_MatMVPY(parms_Mat self, FLOAT alpha, FLOAT *x, FLOAT
			 beta, FLOAT *y, FLOAT *z); 

/** 
 * Perform the multiplication of the off-diagonal matrix and the
 * external vars. 
 *
 * The local matrix can be written as follows:
 * 
 *  \f[
 *  \left(
 *  \begin{array}{ccc}
 *    B   &   E & 0\\
 *    F   &   C & M_{ext}
 *  \end{array}
 *  \right)
 *  \f],
 *  where \f$\left(\begin{array}{cc}
 *    B  &   E \\
 *    F  &   C 
 *  \end{array}
 *  \right)\f$ corresponds to the variables on the local PE. This
 *  function 
 *  performs 
 *  \f[
 *    y[pos..n] = M_{ext} \times x_{ext}
 *  \f]
 *
 * @param self A matrix object.                                      
 * @param x    A vector object.                                      
 * @param y    A vector object.                                      
 * @param pos  The offset of x from the beginning of he local vector.
 * 
 * @return 0 on success.
 */
extern int parms_MatVecOffDiag(parms_Mat self, FLOAT *x, FLOAT *y, int
			       pos);

/** 
 * Get the communication handler.
 * 
 * @param self    A matrix object.                     
 * @param handler The communication handler returned.  
 * 
 * @return 0 on success.
 */
extern int parms_MatGetCommHandler(parms_Mat self, parms_Comm
				   *handler);

/** 
 * Free the memory for the submatrix.
 * 
 * @param self A parms_Mat object.             
 * @param mat  The submatrix to be freed.      
 * 
 * @return 0 on success.
 */
extern int parms_MatFreeSubMat(parms_Mat self, void *mat);

/** 
 * Get the local matrix. 
 * 
 * @param self A matrix object.                               
 * @param mat  The submatrix returned in a specific format.   
 * 
 * @return 0 on success.
 */
extern int parms_MatGetSubMat(parms_Mat self, void **mat);

/** 
 * Extend submatrix by including equations correspond to the
 * immediate neighbouring variables.
 * 
 * @param self     A matrix object.                                         
 * @param handler  A communication handler.                                 
 * @param start    The beginning location of mat in the local matrix.       
 * @param mat      The submatrix to be extended.        
 * @param n 	   The size of extended matrix returned.
 * @param ext_mat  The extended matrix created.            
 * 
 * @return 0 on success.
 */
extern int parms_MatExtend(parms_Mat self, parms_Comm handler, int
			   start, void *mat, int *n, void **ext_mat);

/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_matvec_(parms_Mat *self, FLOAT *x, FLOAT *y, int *ierr);  

extern void parms_matcreate_(parms_Mat *self, parms_Map *map, int *ierr);

extern void parms_matfree_(parms_Mat *self, int *ierr); 

extern void parms_matmvpy_(parms_Mat *self, FLOAT *alpha, FLOAT *x, FLOAT
		    *beta, FLOAT *y, FLOAT *z, int *ierr);

extern void parms_matsetcommtype_(parms_Mat *self, COMMTYPE *ctype, int
			   *ierr);

extern void parms_matsetvalues_(parms_Mat *self, int *m, int *im, int *ia,
			 int *ja, FLOAT *values, INSERTMODE *mode, int
			 *ierr);
			 
extern void parms_matsetelementmatrix_(parms_Mat *self, int *m, int *im, int *ia,
			 int *ja, FLOAT *values, INSERTMODE *mode, int *ierr);

extern void parms_matassembleelementmatrix_(parms_Mat *self, int *ierr);	

extern void parms_matresetrowvalues_(parms_Mat *self, int *m, int *im, int *ia,
		       int *ja, FLOAT *values, int *ierr);	       

extern void parms_matreset_(parms_Mat *self, NNZSTRUCT *nonzerostructure, int *ierr);

extern void parms_matsetup_(parms_Mat *self, int *ierr);

extern void parms_matview_(parms_Mat *self, parms_Viewer *v, int *ierr);

extern void parms_matviewcoo_(parms_Mat *self, parms_Viewer *v, int *ierr);         

/*
 *
 * end Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif 
