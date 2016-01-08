/*!
  \file   parms_mat_impl.h
  \brief  The matrix structure and functions

  \author zzli
  \date   2006-05-08
*/

#ifndef _PARMS_MAT_IMPL_H_
#define _PARMS_MAT_IMPL_H_

#include "parms_mat.h"
#include "parms_operator.h"
#include "parms_map_impl.h"
#include "parms_table_impl.h"

#define PARMS_MAT_FORMAT_KIND    0x00f00000
#define PARMS_MAT_SHIFT          20
#define PARMS_MAT_GET_KIND(a)       (((a)&PARMS_MAT_FORMAT_KIND) >> PARMS_MAT_SHIFT)
#define PARMS_MAT_SET_KIND(a, kind) ((a)|(kind) << PARMS_MAT_SHIFT)

/*
typedef  int (*PARMS_ILU)(parms_Mat self, parms_FactParam param, void
			  *data, parms_Operator *op);  
*/

/*! \struct parms_Mat_ops.
  \brief struct parms_mat_ops.
 */
typedef struct parms_Mat_ops {
  /*! function parameter: apply
    \brief performs \f$y = self \times x\f$. 
  */
  int (*apply)(parms_Mat self, FLOAT *x, FLOAT *y);
  /*! 
    Divide the local equation into diagonal part and off-diagonal
    part. Set up communcation handler for matrix-vector product.
   */
  int (*setup)(parms_Mat self);	
  /*! Set the communication type
    \f{tabular}{ll}
    P2P & Copying data into a temporary buffer \\
    DERIVED & Creating a derived datatype as in the old version of
    pARMS. 
    \f}
   */
  int (*setcommtype)(parms_Mat self, COMMTYPE ctype);

  /*! function parameter: mvpy
    \brief Performs z = beta*y + alpha*self*x.
  */
  int (*mvpy)(parms_Mat self, FLOAT alpha, FLOAT *x, FLOAT beta,
	      FLOAT *y, FLOAT *z);

  /*! function parameter: getdiag
    \brief Get the diagonal part of the local matrix.
   */
  int (*getdiag)(parms_Mat self, void **mat);

  
  /*! function parameter: getoffdiag (used by pc_schurras) (AF)
    \brief Get the offdiagonal part of the local matrix.
   */
  int (*getoffdiag)(parms_Mat self, void **mat);
  /*! function parameter: getlmat
    \brief Get the local matrix including diagonal and off-diagonal
    matrix.  
   */
 
  int (*getlmat)(parms_Mat self, void **mat);
  /*! function parameter: extend
    \brief Extend the submatrix mat by including equations that
    correspond to immediate neighbouring variables.
   */
  int (*extend)(parms_Mat self, parms_Comm handler, int start, void
		*mat, int *n, void **ext_mat);
  /*! function Parameter: mvoffd
    \brief The matrix-vector product for off-diagonal part.
   */
  int (*mvoffd)(parms_Mat self, FLOAT *x, FLOAT *y, int pos);

  /*! function Parameter: matfree
    \brief Free the matrix mat.
   */
  int (*matfree)(parms_Mat self, void *mat);

  /*! function Parameter: gethandler
      \brief Get the communication handler for mv product.
   */
  int (*gethandler)(parms_Mat self, parms_Comm *handler);
} *parms_Mat_ops;


/*! \struct parms_vcsr.
  \ Description: struct parms_vcsr.
 */
typedef struct parms_vcsr {
  int      n;         	//!< the dimension of the matrix 
  int      *nnzrow;    //!< the length of each row 
  int      *space;     //!<  length of space ( a work array)
  int      off_proc_n; //!< size of off-processor contributions
  /*! 
  		Parameter: pj = An indirect pointer to store column indices.
   */
  int      **pj;
  /*! parameter: pa = An indirect pointer to store corresponding nonzero entries.  
  */
  FLOAT    **pa;

  int      *pj_data;
  FLOAT    *pa_data;
} *parms_vcsr;

/*! \struct parms_Mat_
  \Description: struct parms_Mat_
 */
struct parms_Mat_ {

  int ref;
  parms_Mat_ops ops;
  void        *data;
  BOOL        isserial;
  BOOL        issetup;
  BOOL        isperm;
  BOOL        isalloc;
  BOOL        isreset;
  NNZSTRUCT   resetpattern;
  MATTYPE     type;
  PCILUTYPE   ilutype;
  int         m,n;
  int         M,N;
  parms_Map   is;
  parms_vcsr  aux_data;
/* data for external contributions */
  parms_vcsr ext_data;
  BOOL       isassembled;
/* Table to track count of offdiagonal variables.
 * This allows efficient update of the hash table 
 * for the variables, when a row is reset.
*/  
  parms_Table odtable;
};


/* external function protos */
extern int parms_MatCreate_vcsr(parms_Mat self);
extern int parms_MatCreate_dvcsr(parms_Mat self);
extern int parms_MatFree_dvcsr(parms_Mat *self);
extern int parms_MatView_vcsr(parms_Mat self, parms_Viewer v);
extern int parms_MatViewCOO_vcsr(parms_Mat self, parms_Viewer v);
extern int parms_MatView_dvcsr(parms_Mat self, parms_Viewer v);
extern int parms_MatViewCOO_dvcsr(parms_Mat self, parms_Viewer v);

#endif 
