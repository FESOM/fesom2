/**
 * @file   parms_vec.h
 * @author Zhongze Li
 * @date   Tue Oct 17 12:01:25 2006
 * 
 * @brief  Functions related to the parms_Vec object. 
 * 
 * 
 */

#ifndef _PARMS_VECTOR_H_
#define _PARMS_VECTOR_H_

#include "parms_sys.h"
#include "parms_map.h"
#include "parms_viewer.h"
#include "parms_comm.h"

PARMS_CXX_BEGIN

/** 
 * Return the 2-norm of the vector.
 * 
 * @param self  A vector object.    
 * @param value The 2-norm returned.
 * 
 * @return 0 on success.
 */
extern int parms_VecGetNorm2(FLOAT *self, REAL *value, parms_Map map);

/** 
 * Scale a vector.
 *
 * All components of vector self on the local processor times scalar. 
 * \f$self = scalar \times self\f$.
 *  
 * @param self   A vector object.
 * @param scalar A scalar.      
 * 
 * @return 0 on success.
 */
extern int parms_VecScale(FLOAT *self, FLOAT scalar, parms_Map map);

/** 
 * Perform \f$self := scalar \times x + self\f$.
 * 
 * @param self   A vector object.      
 * @param x      Another vector object.
 * @param scalar A scalar.
 * 
 * @return 0 on success.
 */
extern int parms_VecAXPY(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map);

/** 
 * Perform \f$self = scalar \times self + x\f$.
 * 
 * @param self    A vector object.      
 * @param x       Another vector object.
 * @param scalar  A scalar.
 * 
 * @return 0 on success.
 */
extern int parms_VecAYPX(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map);

/** 
 * Perform the (global) inner product of two vectors.
 * 
 *  value = self x^{T}. If self
 * 
 *
 * @param self   A vector object.                  
 * @param x 	 Another vector object.            
 * @param value  The inner product returned.       
 * 
 * @return 0 on success.
 */
extern int parms_VecDOT(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map);

/** 
 * Perform the (global) inner product of two vectors.
 * 
 * If self and x are real vectors, value = self x^{T}. If self
 * and x are complex vectors, value = self \overline{x}^{T}.
 *
 * @param self   A vector object.                  
 * @param x 	 Another vector object.            
 * @param value  The inner product returned.       
 * 
 * @return 0 on success.
 */
int parms_VecDOTC(FLOAT *self, FLOAT *x, REAL *value, parms_Map is);

/** 
 * Perform the inner product between self and an array of parms_Vec
 * objects.  
 *
 * The pseudo code:
 *
 *  \f{verbatim}
 *  for (i = 0; i < n; i++) {
 *    result[i] = self * vecarray[i];
 *  }
 *  \f}
 *
 * @param self     A vector object.                              
 * @param n        The size of vecarray.                         
 * @param vecarray An array of vector objects.                   
 * @param aux      An auxiliary array.                           
 * @param result   An array of size n to store inner products.   
 * 
 * @return 0 on success.
 */
extern int parms_VecDotArray(FLOAT *self, int n, FLOAT
			     **vecarray, FLOAT *result, parms_Map map);   
/** 
 * Permute the vector object self.
 *
 * If the parms_Vec object and the parms_Mat object are created based
 * on  the same parms_Map object. Once matrix object is setup,
 * variables on each processor are divided into two cateries: internal
 * unknowns and interface unknowns. The parms_Vec should be permuted
 * accordingly. The user needn't call self function directly.
 * 
 * @param self A vector object.
 * 
 * @return 0 on success.
 */
extern int parms_VecPerm(FLOAT *self, parms_Map map);

/** 
 * Inverse permutation of the vector object self.
 * The user needn't call this function directly.
 * 
 * @param self A vector object.
 * 
 * @return 0 on success.
 */
extern int parms_VecInvPerm(FLOAT *self, parms_Map map);

/** 
 * Permute the vector object self into the vector aux.
 *
 * This uses the local permutation based on the local map object.
 * The user needn't call self function directly.
 * 
 * @param self A vector object.
 * @param aux  A vector object containing the permuted vector on return.
 * 
 * @return 0 on success.
 */
extern int parms_VecPermAux(FLOAT *self, FLOAT *aux, parms_Map map);

/** 
 * Inverse permutation the vector object self into the vector aux.
 *
 * This uses the local permutation based on the local map object.
 * The user needn't call self function directly.
 * 
 * @param self A vector object.
 * @param aux  A vector object containing the (inverse) permuted vector on return.
 * 
 * @return 0 on success.
 */
extern int parms_VecInvPermAux(FLOAT *self, FLOAT *aux, parms_Map map);

/** 
 * Insert or add values to vector object self.
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables to be inserted.  
 * @param im     An array of global variable indices.        
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
extern int parms_VecSetValues(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map); 

/** 
 * Insert values to parms_Vec object self. This assumes the vector 
 * values are being set element-by-element. A call to parms_vecsetupElements
 * is required to complete the vector once all entries have been added.
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables to be inserted.     
 * @param im     An array of global variable indices.     
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
int parms_VecSetElementVector(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map);
			       
/** 
 * Completes setting up values for the distributed vector
 *
 * @param self   A vector object.                          
 * @param map    A pARMS map object
 *               
 * @return 0 on success.
 */
int parms_VecAssembleElementVector(FLOAT *self, parms_Map map); 			        

/** 
 * Gather distributed vector to a global array.
 * 
 * @param self The distributed vector.
 * @param ga A global vector.
 * 
 * @return 0 on success.
 */
extern int parms_VecGather(FLOAT *self, FLOAT *ga, parms_Map map);

/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_vecaxpy_(FLOAT *self, FLOAT *x, FLOAT *scalar, parms_Map *map, int
		    *ierr);

extern void parms_vecaypx_(FLOAT *self, FLOAT *x, FLOAT *scalar, parms_Map *map, int
		    *ierr);

extern void parms_vecdot_(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map *map, int
		   *ierr);

extern void parms_vecdotc_(FLOAT *self, FLOAT *x, REAL *value, parms_Map *map, int
		   *ierr);

extern void parms_vecdotarray_(FLOAT *self, int *n, FLOAT **vecarray,
			FLOAT *result, parms_Map *map, int *ierr);

extern void parms_vecgetnorm2_(FLOAT *self, REAL *value, parms_Map *map, int *ierr);

extern void parms_vecscale_(FLOAT *self, FLOAT *scalar, parms_Map *map, int *ierr);

extern void parms_vecsetvalues_(FLOAT *self, int *m, int *im, FLOAT *values, 
                  INSERTMODE *mode, parms_Map *map, int *ierr);
                  
extern void parms_vecsetelementvector_(FLOAT *self, int *m, int *im, FLOAT *values, 
                  INSERTMODE *mode, parms_Map *map, int *ierr); 
			       
extern void parms_vecassembleelementvector_(FLOAT *self, parms_Map *map, int *ierr);

extern void parms_vecperm_(FLOAT *self, parms_Map *map, int *ierr);

extern void parms_vecinvperm_(FLOAT *self, parms_Map *map, int *ierr);

extern void parms_vecpermaux_(FLOAT *self, FLOAT *aux, parms_Map *map, int *ierr);

extern void parms_vecinvpermaux_(FLOAT *self, FLOAT *aux, parms_Map *map, int *ierr);

extern void parms_vecgather_(FLOAT *self, FLOAT *ga, parms_Map *map, int *ierr);


/*
 *
 * End Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif 
