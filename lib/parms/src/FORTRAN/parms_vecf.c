#include "parms_vec.h"

#if defined(FORTRAN_CAPS)
#define parms_vecaxpy_                PARMS_VECAXPY            
#define parms_vecaypx_		      PARMS_VECAYPX            
#define parms_vecdot_		      PARMS_VECDOT             
#define parms_vecdotc_		      PARMS_VECDOTC             
#define parms_vecdotarray_	      PARMS_VECDOTARRAY        
#define parms_vecgetnorm2_	      PARMS_VECGETNORM2        
#define parms_vecscale_		      PARMS_VECSCALE           
#define parms_vecsetvalues_	      PARMS_VECSETVALUES       
#define parms_vecsetelementvector_	      PARMS_VECSETELEMENTVECTOR       
#define parms_vecassembleelementvector_	      PARMS_VECASSEMBLEELEMENTVECTOR       
#define parms_vecperm_             PARMS_VECPERM
#define parms_vecinvperm_             PARMS_VECINVPERM
#define parms_vecinvpermaux_          PARMS_VECINVPERMAUX
#define parms_vecpermaux_          PARMS_VECPERMAUX
#define parms_vecgather_              PARMS_VECGATHER
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_vecaxpy_                parms_vecaxpy__            
#define parms_vecaypx_		      parms_vecaypx__            
#define parms_vecdot_		      parms_vecdot__             
#define parms_vecdotc_		      parms_vecdotc__             
#define parms_vecdotarray_	      parms_vecdotarray__        
#define parms_vecgetnorm2_	      parms_vecgetnorm2__        
#define parms_vecscale_		      parms_vecscale__           
#define parms_vecsetvalues_	      parms_vecsetvalues__       
#define parms_vecsetelementvector_	      parms_vecsetelementvector__       
#define parms_vecassembleelementvector_	      parms_vecassembleelementvector__       
#define parms_vecperm_             parms_vecperm__
#define parms_vecinvperm_             parms_vecinvperm__
#define parms_vecinvpermaux_          parms_vecinvpermaux__
#define parms_vecpermaux_          parms_vecpermaux__
#define parms_vecgather_              parms_vecgather__
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_vecaxpy_                parms_vecaxpy            
#define parms_vecaypx_		      parms_vecaypx            
#define parms_vecdot_		      parms_vecdot             
#define parms_vecdotc_		      parms_vecdotc             
#define parms_vecdotarray_	      parms_vecdotarray        
#define parms_vecgetnorm2_	      parms_vecgetnorm2        
#define parms_vecscale_		      parms_vecscale           
#define parms_vecsetvalues_	      parms_vecsetvalues       
#define parms_vecsetelementvector_	      parms_vecsetelementvector       
#define parms_vecsassembleelementvector_	      parms_vecassembleelementvector       
#define parms_vecperm_             parms_vecperm
#define parms_vecinvperm_             parms_vecinvperm
#define parms_vecinvpermaux_          parms_vecinvpermaux
#define parms_vecpermaux_          parms_vecpermaux
#define parms_vecgather_              parms_vecgather
#endif

void parms_vecaxpy_(FLOAT *self, FLOAT *x, FLOAT *scalar, parms_Map *map, int
		    *ierr)
{
  *ierr = parms_VecAXPY(self, x, *scalar, *map);
}

void parms_vecaypx_(FLOAT *self, FLOAT *x, FLOAT *scalar, parms_Map *map, int
		    *ierr)
{
  *ierr = parms_VecAYPX(self, x, *scalar, *map);
}


void parms_vecdot_(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map *map, int
		   *ierr) 
{
  *ierr = parms_VecDOT(self, x,  value, *map);
}

void parms_vecdotc_(FLOAT *self, FLOAT *x, REAL *value, parms_Map *map, int
		   *ierr) 
{
  *ierr = parms_VecDOTC(self, x,  value, *map);
}

void parms_vecdotarray_(FLOAT *self, int *n, FLOAT **vecarray,
			FLOAT *result, parms_Map *map, int *ierr) 
{
  *ierr = parms_VecDotArray(self, *n, vecarray, result, *map);   
}

void parms_vecgetnorm2_(FLOAT *self, REAL *value, parms_Map *map, int *ierr)
{
  *ierr = parms_VecGetNorm2(self, value, *map);
}

void parms_vecscale_(FLOAT *self, FLOAT *scalar, parms_Map *map, int *ierr)
{
  *ierr = parms_VecScale(self, *scalar, *map);
}

void parms_vecsetvalues_(FLOAT *self, int *m, int *im, FLOAT *values, INSERTMODE *mode, parms_Map *map, int *ierr)
{
  *ierr = parms_VecSetValues(self, *m, im, values, *mode, *map) ;
}

void parms_vecsetelementvector_(FLOAT *self, int *m, int *im, FLOAT *values, INSERTMODE *mode, parms_Map *map, int *ierr)
{
  *ierr = parms_VecSetElementVector(self, *m, im, values, *mode, *map) ;
}

void parms_vecassembleelementvector_(FLOAT *self, parms_Map *map, int *ierr)
{
  *ierr = parms_VecAssembleElementVector(self, *map) ;
}

void parms_vecperm_(FLOAT *self, parms_Map *map, int *ierr)
{
  *ierr = parms_VecPerm(self, *map);
}

void parms_vecinvperm_(FLOAT *self, parms_Map *map, int *ierr)
{
  *ierr = parms_VecInvPerm(self, *map);
}

void parms_vecpermaux_(FLOAT *self, FLOAT *aux, parms_Map *map, int *ierr)
{
  *ierr = parms_VecPermAux(self, aux, *map);
}

void parms_vecinvpermaux_(FLOAT *self, FLOAT *aux, parms_Map *map, int *ierr)
{
  *ierr = parms_VecInvPermAux(self, aux, *map);
}

void parms_vecgather_(FLOAT *self, FLOAT *ga, parms_Map *map, int *ierr)
{
  *ierr = parms_VecGather(self, ga, *map);
}


