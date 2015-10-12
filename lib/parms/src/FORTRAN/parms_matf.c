#include "parms_mat.h"

#if defined(FORTRAN_CAPS)
#define parms_matvec_           PARMS_MATVEC        
#define parms_matcreate_ 	  PARMS_MATCREATE       
#define parms_matfree_ 		  PARMS_MATFREE         
#define parms_matmvpy_ 		  PARMS_MATMVPY         
#define parms_matsetcommtype_ 	  PARMS_MATSETCOMMTYPE    
#define parms_matsetvalues_ 	  PARMS_MATSETVALUES 
#define parms_matsetelementmatrix_ 	  PARMS_MATSETELEMENTMATRIX    
#define parms_matassembleelementmatrix_ 	  PARMS_MATASSEMBLEELEMENTMATRIX    
#define parms_matresetrowvalues_ 	  PARMS_MATRESETROWVALUES    
#define parms_matreset_ 	  PARMS_MATRESET    
#define parms_matsetup_ 	  PARMS_MATSETUP        
#define parms_matview_         	  PARMS_MATVIEW         
#define parms_matviewcoo_         	  PARMS_MATVIEWCOO         
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_matvec_           parms_matvec__        
#define parms_matcreate_ 	  parms_matcreate__       
#define parms_matfree_ 		  parms_matfree__         
#define parms_matmvpy_ 		  parms_matmvpy__         
#define parms_matsetcommtype_ 	  parms_matsetcommtype__     
#define parms_matsetvalues_ 	  parms_matsetvalues__    
#define parms_matsetelementmatrix_ 	  parms_matsetelementmatrix__    
#define parms_matassembleelementmatrix_ 	  parms_matassembleelementmatrix__    
#define parms_matresetrowvalues_ 	  parms_matresetrowvalues__    
#define parms_matreset_ 	  parms_matreset__    
#define parms_matsetup_ 	  parms_matsetup__        
#define parms_matview_         	  parms_matview__         
#define parms_matviewcoo_         	  parms_matviewcoo__         
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_matvec_           parms_matvec        
#define parms_matcreate_ 	  parms_matcreate       
#define parms_matfree_ 		  parms_matfree         
#define parms_matmvpy_ 		  parms_matmvpy         
#define parms_matsetcommtype_ 	  parms_matsetcommtype       
#define parms_matsetvalues_ 	  parms_matsetvalues    
#define parms_matsetelementmatrix_ 	  parms_matsetelementmatrix    
#define parms_matassembleelementmatrix_ 	  parms_matassembleelementmatrix    
#define parms_matresetrowvalues_ 	  parms_matresetrowvalues    
#define parms_matreset_ 	  parms_matreset    
#define parms_matsetup_ 	  parms_matsetup        
#define parms_matview_         	  parms_matview         
#define parms_matviewcoo_         	  parms_matviewcoo         
#endif

void parms_matvec_(parms_Mat *self, FLOAT *x, FLOAT *y, int
		     *ierr)  
{
  *ierr = parms_MatVec(*self, x, y);
}

void parms_matcreate_(parms_Mat *self, parms_Map *map, int *ierr)        
{
  *ierr = parms_MatCreate(self, *map);
}

void parms_matfree_(parms_Mat *self, int *ierr)          
{
  *ierr = parms_MatFree(self);
}

void parms_matmvpy_(parms_Mat *self, FLOAT *alpha, FLOAT *x, FLOAT
		    *beta, FLOAT *y, FLOAT *z, int *ierr)
{
  *ierr = parms_MatMVPY(*self, *alpha, x, *beta, y, z);
}

void parms_matsetcommtype_(parms_Mat *self, COMMTYPE *ctype, int
			   *ierr)
{
  *ierr = parms_MatSetCommType(*self, *ctype);
}

void parms_matsetvalues_(parms_Mat *self, int *m, int *im, int *ia,
			 int *ja, FLOAT *values, INSERTMODE *mode, int
			 *ierr)  
{
  *ierr = parms_MatSetValues(*self, *m, im, ia, ja, values, *mode);
}

void parms_matsetelementmatrix_(parms_Mat *self, int *m, int *im, int *ia,
			 int *ja, FLOAT *values, INSERTMODE *mode, int
			 *ierr)  
{
  *ierr = parms_MatSetElementMatrix(*self, *m, im, ia, ja, values, *mode);
}

void parms_matassembleelementmatrix_(parms_Mat *self, int *ierr)  
{
  *ierr = parms_MatAssembleElementMatrix(*self);
}

void parms_matresetrowvalues_(parms_Mat *self, int *m, int *im, int *ia,
		       int *ja, FLOAT *values, int *ierr)
{
  *ierr = parms_MatResetRowValues(*self, *m, im, ia, ja, values);
}

void parms_matreset_(parms_Mat *self, NNZSTRUCT *nonzerostructure, int *ierr)		       
{
  *ierr = parms_MatReset(*self, *nonzerostructure);
}

void parms_matsetup_(parms_Mat *self, int *ierr)
{
  *ierr = parms_MatSetup(*self);
}

void parms_matview_(parms_Mat *self, parms_Viewer *v, int *ierr)        
{
  *ierr = parms_MatView(*self, *v);
}

void parms_matviewcoo_(parms_Mat *self, parms_Viewer *v, int *ierr)        
{
  *ierr = parms_MatViewCOO(*self, *v);
}

