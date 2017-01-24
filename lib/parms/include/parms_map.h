/**
 * @file   parms_map.h
 * @author Zhongze Li
 * @date   Tue Oct 17 10:02:53 2006
 * 
 * @brief  Functions related to the parms_Map object. 
 * 
 * parms_Map describes how variables are distributed across
 * processors. It is used for creating parms_Mat
 * objects. 
 */

#ifndef _PARMS_MAP_H_
#define _PARMS_MAP_H_

#include "parms_sys.h"
#include "parms_viewer.h"

PARMS_CXX_BEGIN

typedef struct parms_Map_ *parms_Map;

/** 
 * Create a parms_Map object on the local processor.
 * 
 * @param self   A pointer to the parms_Map object created.
 * @param gsize  The size of unknowns on the local processor.       
 * @param offset The start index.
 *               - 1 FORTRAN
 *               - 0 C
 * 
 * @return 0 on success.
 *
 */
extern int parms_MapCreateFromLocal(parms_Map *self, int gsize, int
				    offset);  

/** 
 * Create a parms_Map object based on the Metis partitioning.
 * 
 * @param self     A pointer to the parms_Map object created.       
 * @param gsize    The total number of vertices.                    
 * @param npar 	   An integer array of size gsize. node \f$i\f$     
 *       	   resides on PE npar[i].                           
 * @param comm     MPI communicator.                                
 * @param offset   The start index.                                 
 *                 - 1 FORTRAN
 *                 - 0 C
 * @param dof      The number of variables associated with each
 *                 vertex.   
 * @param VARSTYPE Assuming the variables \f$u_i, v_i\f$ are
 *                 associated with vertex \f$i\f$, two styles of
 *                 numbering variables are as follows:
 *                 - INTERLACED. Variables are numbered in the
 *                   order of \f$u_1, v_1, u_2, v_2, \cdots\f$; 
 *                 - NONINTERLACED. Variables are numbered in the
 *                   order of \f$u_1, u_2, u_3,...,v_1, v_2,...\f$. 
 *
 * @return 0 on success.
 */
extern int parms_MapCreateFromGlobal(parms_Map *self, int gsize, int *npar,  
				     MPI_Comm comm, int offset, int dof, 
				     VARSTYPE vtype);  
				     
/** 
 * Create a parms_Map object based on the output of ParMetis.
 * 
 * @param self     A parms_Map object created.
 * @param vtxdist  An integer array of size np+1, where np is the
 *                 number of PEs. This array indicates the range of
 *                 vertices that are local to each processor.  PE i
 *                 stores vertices in the range of [vtxdist[i],
 *                 vtxdist[i+1]).
 * @param part     An array of size equal to the number of
 *                 locally-stored vertices. part[j] indicates the ID
 *                 of the PE to which the vertex with local index j
 *                 and global index vtxdist[pid]+j belongs (pid is ID
 *                 of local PE). 
 * @param comm     MPI communicator.
 * @param offset   The start index. 
 *                 - 1 FORTRAN
 *                 - 0 C
 * @param dof      The number of variables associated with each
 *                 vertex.  
 * @param vtype    Assuming the variables u_i, v_i are associated
 *                 with vertex i, two styles of numbering variables
 *                 are as follows: 
 *                 - INTERLACED. Variables are numbered in the
 *                   order of \f$u_1, v_1, u_2, v_2, \cdots\f$; 
 *                 - NONINTERLACED. Variables are numbered in the
 *                   order of \f$u_1, u_2, u_3,...,v_1, v_2,...\f$. 
 * 
 * @return 0 on success.
 */
extern int parms_MapCreateFromDist(parms_Map *self, int *vtxdist, int
				   *part, MPI_Comm comm, int offset, 
				   int dof, VARSTYPE vtype);   

/** 
 * Create a parms_Map object with the default partioning strategy in
 * PETSc. 
 * 
 * @param self  A parms_Map object created.  
 * @param m     The local size of variables. 
 * @param M     The global size of variables.
 * @param comm  MPI communicatior.           
 * 
 * @return 0 on success.
 */
extern int parms_MapCreateFromPetsc(parms_Map *self, int m, int M,
				    MPI_Comm comm);  

/** 
 * Create a parms_Map object.
 * 
 * @param self     A parms_Map object created.
 * @param gsize    The total number of vertices.       
 * @param nodes    A list of all vertices stored PE by PE.        
 * @param p2nodes  An integer array of size np+1, np is the number of
 *                 PEs. If k1 = p2nodes[i], k2 = p2nodes[i+1]
 *                 , then PE i contains the vertices in the 
 *                 range of [[nodes[k1], [nodes[k2-1]]. 
 * @param comm     MPI communication.
 * @param dof      The number of variables associated with each
 *                 node. 
 * @param vtype    Assuming the variables \f$u_i, v_i\f$ are
 *                 associated with vertex \f$i\f$, two style of
 *                 numbering variables are as follows:
 *                 - INTERLACED. Variables are numbered in the
 *                   order of \f$u_1, v_1, u_2, v_2, \cdots\f$; 
 *                 - NONINTERLACED. Variables are numbered in the
 *                   order of \f$u_1, u_2, u_3,...,v_1, v_2,...\f$. 
 * 
 * @return 0 on success.
 */
extern int parms_MapCreateFromPtr(parms_Map *self, int gsize, int
				  *nodes, int *p2nodes, MPI_Comm 
				  comm, int dof, VARSTYPE vtype);  

/** 
 * Free the memory for the parms_Map object pointed to by self.
 * 
 * @param self A pointer to the parms_Map object.
 * 
 * @return 0 on success.
 */
extern int parms_MapFree(parms_Map *self);

/** 
 * Get the size of variables on the local PE.
 * 
 * @param self A parms_Map object.
 * 
 * @return The size of variables on the local PE.
 */
extern int parms_MapGetLocalSize(parms_Map self);

/** 
 * Get global size of variables.
 * 
 * @param self A parms_Map object.
 * 
 * @return  The global size of variables rather than vertices.
 */
extern int parms_MapGetGlobalSize(parms_Map self);

/** 
 * Get PE's ID.
 * 
 * @param self A parms_Map object.
 * 
 * @return The PE's ID.
 */
extern int parms_MapGetPid(parms_Map self);

/** 
 * Get the number of PEs.
 * 
 * @param self A parms_Map object.
 * 
 * @return The number of PEs.
 */
extern int parms_MapGetNumProcs(parms_Map self);

/** 
 * Dump the parms_Map object.
 * 
 * @param self A parms_Map object.
 * @param v    A parms_Viewer object.    
 * 
 * @return 0 on success.
 */
extern int parms_MapView(parms_Map self, parms_Viewer v);

/** 
 * Get local index for a given global index.
 *
 * Return a pointer to an integer. If it is NULL, then the variable
 * with global index gindex doesn't reside on the local
 * PE. Otherwise, it points to an address of a variable 
 * whose value is local index. 
 *
 * @param self    A parms_Map object.
 * @param gindex  A global index.
 * 
 * @return A pointer to an integer whose value is the corresponding
 *         local index. 
 */
extern int *parms_MapGlobalToLocal(parms_Map self, int gindex);

/**
 * Get local variable array
 *
*/
extern int parms_MapGetGlobalIndices(parms_Map self, int *im);

/**
 *
 * Fortran Wrapper Functions
 *
 *
*/
extern void parms_mapcreatefromlocal_(parms_Map *self, int *gsize, int
			       *offset, int *ierr);  

extern void parms_mapcreatefromglobal_(parms_Map *self, int *gsize, int
				*npar, MPI_Comm *comm, int *offset,
				int *dof, VARSTYPE *vtype, int *ierr);  
				
extern void parms_mapcreatefromdist_(parms_Map *self, int *vtxdist, int
			     *part, MPI_Comm *comm, int *offset, int
			     *dof, VARSTYPE *vtype, int *ierr);
			     
extern void parms_mapcreatefrompetsc_(parms_Map *self, int *m, int *M,
			       MPI_Comm *comm, int *ierr); 

extern void parms_mapcreatefromptr_(parms_Map *self, int *gsize, int *nodes,
			     int *p2nodes, MPI_Comm *comm, int *dof,
			     VARSTYPE *vtype, int *ierr); 

extern void parms_mapfree_(parms_Map *self, int *ierr);

extern void parms_mapgetglobalsize_(parms_Map *self, int *gsize);

extern void parms_mapgetlocalsize_(parms_Map *self, int *lsize);

extern void parms_mapgetnumprocs_(parms_Map *self, int *numpro);

extern void parms_mapgetpid_(parms_Map *self, int *pid);

extern void parms_mapview_(parms_Map *self, parms_Viewer *v, int *ierr);

extern void parms_mapgetglobalindices_(parms_Map *self, int *im, int *ierr);
/*
 *
 * end Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif 
