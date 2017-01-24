/*!
  \file   parms_map_impl.h
  \brief  The definition of the parms_Map struct.

  \author zzli
  \date   2006-05-05
*/

#ifndef _PARMS_MAP_IMPL_H_
#define _PARMS_MAP_IMPL_H_

#include "parms_mem.h"
#include "parms_table.h"
#include "parms_map.h"


/*! \struct parms_Map_
 */
struct parms_Map_ {
  int ref;
  parms_Table table;		/**< contains pair (gindex,lindex) */
  MPI_Comm    comm;		/**< MPI communicator */
  int         pid;		/**< processor ID */
  int         npro;		/**< number of processors */
  int         lsize;		/**< number of local variables */
  int         gsize;		/**< number of global variables */
  /*! numbering style used.
   *
   *  \f{tabular}{lcl} 
   *   FORTRAN &-& 1 \\
   *   C       &-& 0 
   *   \f}
   */
  int         start;
  /*! 
     The number of variables associated with each vertex.
  */
  int         dof;	
  /*! style of labelling variables u_i,v_i associated with the i-th vertex.
   *
   * - NONINTERLACED u_1,v_1,u_2,v_2.
   *
   * - INTERLACED u_1,u_2,\cdots,v_1,v_2.
   */
  VARSTYPE    vtype;	
  /*! array of size lsize, stores local variables in global labels */
  int         *lvars;		
  BOOL        isserial;		//!< data layout on one processor?
  /*! variables are permuted or not.
   *  -true: interior variables labelled first followed by interface
   *  variables */
  BOOL        isperm;		
  BOOL	      isvecperm;	/* used to check if vector has been permuted in gmres */	
  BOOL        ispermalloc;	/* permutation array allocated? */
  int         *perm;		//!< permutation array of size lsize.
  int         *iperm;		//!< inverse permutation array. 
  int         nint;		//!< number of interior variables.
  int         ninf;		//!< number of interface variables.
  /*! 
    \param schur_start  start of the local Schur complement which
    may be lesser than nbnd for the matrix with unsymmetric pattern 
  */
  int         schur_start;	
  int         ninf_send;	//!< number of variables to be sent.
  /*! 
    \param vsend  array of size ninf_send which stores variables in local
         indices 
  */
  int         *vsend;	
  /*! 
    \parameters for external data contributions
  */
  int n_ext;  
  int *ext_im; 
  /*! 
    \parameter data used to handle external vector contributions
  */     
  void *data;
  BOOL isdatalloc;
  /*! 
    \param vstable  stores variables to be sent to adjacent processors in pairs
   	 (local_index, index_in_vsend) 
   */
  parms_Table vstable;		
};

#endif 
