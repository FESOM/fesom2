#include "parms_map.h"

#if defined(FORTRAN_CAPS)
#define parms_mapcreatefromlocal_    PARMS_MAPCREATEFROMLOCAL  
#define parms_mapcreatefromglobal_   PARMS_MAPCREATEFROMGLOBAL  
#define parms_mapcreatefromdist_     PARMS_MAPCREATEFROMDIST
#define parms_mapcreatefrompetsc_    PARMS_MAPCREATEFROMPETSC  
#define parms_mapcreatefromptr_	     PARMS_MAPCREATEFROMPTR    
#define parms_mapfree_		     PARMS_MAPFREE         
#define parms_mapgetglobalsize_	     PARMS_MAPGETGLOBALSIZE
#define parms_mapgetlocalsize_	     PARMS_MAPGETLOCALSIZE 
#define parms_mapgetnumprocs_	     PARMS_MAPGETNUMPROCS  
#define parms_mapgetpid_	     PARMS_MAPGETPID       
#define parms_mapview_         	     PARMS_MAPVIEW       
#define parms_mapgetglobalindices_   PARMS_MAPGETGLOBALINDICES_	
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_mapcreatefromlocal_    parms_mapcreatefromlocal__  
#define parms_mapcreatefromglobal_   parms_mapcreatefromglobal__  
#define parms_mapcreatefromdist_     parms_mapcreatefromdist__
#define parms_mapcreatefrompetsc_    parms_mapcreatefrompetsc__  
#define parms_mapcreatefromptr_	     parms_mapcreatefromptr__    
#define parms_mapfree_		     parms_mapfree __
#define parms_mapgetglobalsize_	     parms_mapgetglobalsize__
#define parms_mapgetlocalsize_	     parms_mapgetlocalsize__
#define parms_mapgetnumprocs_	     parms_mapgetnumprocs__
#define parms_mapgetpid_	     parms_mapgetpid__
#define parms_mapview_         	     parms_mapview__         
#define parms_mapgetglobalindices_   parms_mapgetglobalindices__	
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_mapcreatefromlocal_    parms_mapcreatefromlocal  
#define parms_mapcreatefromglobal_   parms_mapcreatefromglobal  
#define parms_mapcreatefromdist_     parms_mapcreatefromdist
#define parms_mapcreatefrompetsc_    parms_mapcreatefrompetsc  
#define parms_mapcreatefromptr_	     parms_mapcreatefromptr    
#define parms_mapfree_		     parms_mapfree
#define parms_mapgetglobalsize_	     parms_mapgetglobalsize
#define parms_mapgetlocalsize_	     parms_mapgetlocalsize
#define parms_mapgetnumprocs_	     parms_mapgetnumprocs
#define parms_mapgetpid_	     parms_mapgetpid
#define parms_mapview_         	     parms_mapview         
#define parms_mapgetglobalindices_   parms_mapgetglobalindices	
#endif

void parms_mapcreatefromlocal_(parms_Map *self, int *gsize, int
			       *offset, int *ierr)  
{
  *ierr = parms_MapCreateFromLocal(self, *gsize, *offset);
}

void parms_mapcreatefromglobal_(parms_Map *self, int *gsize, int
				*npar, MPI_Comm *comm, int *offset,
				int *dof, VARSTYPE *vtype, int *ierr)  
{
  *ierr = parms_MapCreateFromGlobal(self, *gsize, npar, MPI_COMM_WORLD,
				    *offset, *dof, *vtype);    
}

void parms_mapcreatefromdist_(parms_Map *self, int *vtxdist, int
			     *part, MPI_Comm *comm, int *offset, int
			     *dof, VARSTYPE *vtype, int *ierr)
{
  *ierr = parms_MapCreateFromDist(self, vtxdist, part, MPI_COMM_WORLD, *offset,
			     *dof, *vtype);    
}

void parms_mapcreatefrompetsc_(parms_Map *self, int *m, int *M,
			       MPI_Comm *comm, int *ierr) 
{
  *ierr =  parms_MapCreateFromPetsc(self, *m, *M, MPI_COMM_WORLD); 
}

void parms_mapcreatefromptr_(parms_Map *self, int *gsize, int *nodes,
			     int *p2nodes, MPI_Comm *comm, int *dof,
			     VARSTYPE *vtype, int *ierr) 
{
  *ierr = parms_MapCreateFromPtr(self, *gsize, nodes, p2nodes, MPI_COMM_WORLD,
			     *dof, *vtype);  
}

void parms_mapfree_(parms_Map *self, int *ierr)
{
  *ierr = parms_MapFree(self);
}

void parms_mapgetglobalsize_(parms_Map *self, int *gsize)
{
  *gsize = parms_MapGetGlobalSize(*self);
}

void parms_mapgetlocalsize_(parms_Map *self, int *lsize)
{
  *lsize = parms_MapGetLocalSize(*self);
}

void parms_mapgetnumprocs_(parms_Map *self, int *numpro)
{
  *numpro = parms_MapGetNumProcs(*self);
}

void parms_mapgetpid_(parms_Map *self, int *pid)
{
  *pid = parms_MapGetPid(*self);
}

void parms_mapview_(parms_Map *self, parms_Viewer *v, int *ierr)
{
  *ierr = parms_MapView(*self, *v);
}

void parms_mapgetglobalindices_(parms_Map *self, int *im, int *ierr)
{
  *ierr = parms_MapGetGlobalIndices(*self, im);
}

