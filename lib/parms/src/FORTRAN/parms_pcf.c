#include "parms_mem.h"
#include "parms_pc.h"

#if defined(FORTRAN_CAPS)
#define parms_pccreate_                   PARMS_PCCREATE
#define parms_pcfree_ 			  PARMS_PCFREE
#define parms_pcgetratio_ 		  PARMS_PCGETRATIO
#define parms_pcsetbsize_ 		  PARMS_PCSETBSIZE
#define parms_pcsetfill_ 		  PARMS_PCSETFILL
#define parms_pcsetilutype_ 		  PARMS_PCSETILUTYPE
#define parms_pcsetinnereps_ 		  PARMS_PCSETINNEREPS
#define parms_pcsetinnerksize_ 		  PARMS_PCSETINNERKSIZE
#define parms_pcsetinnermaxits_ 	  PARMS_PCSETINNERMAXITS
#define parms_pcsetnlevels_ 		  PARMS_PCSETNLEVELS
#define parms_pcsetparams_ 		  PARMS_PCSETPARAMS
#define parms_pcsettol_ 		  PARMS_PCSETTOL
#define parms_pcsettolind_ 		  PARMS_PCSETTOLIND
#define parms_pcsettype_ 		  PARMS_PCSETTYPE
#define parms_pcsetup_ 			  PARMS_PCSETUP
#define parms_pcsolve_ 			  PARMS_PCSOLVE
#define parms_pcview_            	  PARMS_PCVIEW
#define parms_pcgetname_                  PARMS_PCGETNAME
#define parms_pcilugetname_               PARMS_PCILUGETNAME
#define parms_pccreateabstract_           PARMS_PCCREATEABSTRACT
#define parms_pcsetop_                    PARMS_PCSETOP
#define parms_pcsetpermscaloptions_       PARMS_PCSETPERMSCALOPTIONS
#define parms_pcsetpermtype_              PARMS_PCSETPERMTYPE
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_pccreate_                   parms_pccreate__
#define parms_pcfree_ 			  parms_pcfree__
#define parms_pcgetratio_ 		  parms_pcgetratio__
#define parms_pcsetbsize_ 		  parms_pcsetbsize__
#define parms_pcsetfill_ 		  parms_pcsetfill__
#define parms_pcsetilutype_ 		  parms_pcsetilutype__
#define parms_pcsetinnereps_ 		  parms_pcsetinnereps__
#define parms_pcsetinnerksize_ 		  parms_pcsetinnerksize__
#define parms_pcsetinnermaxits_ 	  parms_pcsetinnermaxits__
#define parms_pcsetnlevels_ 		  parms_pcsetnlevels__
#define parms_pcsetparams_ 		  parms_pcsetparams__
#define parms_pcsettol_ 		  parms_pcsettol__
#define parms_pcsettolind_ 		  parms_pcsettolind__
#define parms_pcsettype_ 		  parms_pcsettype__
#define parms_pcsetup_ 			  parms_pcsetup__
#define parms_pcsolve_ 			  parms_pcsolve__
#define parms_pcview_            	  parms_pcview__
#define parms_pcgetname_                  parms_pcgetname__
#define parms_pcilugetname_               parms_pcilugetname__
#define parms_pccreateabstract_           parms_pccreateabstract__
#define parms_pcsetop_                    parms_pcsetop__
#define parms_pcsetpermscaloptions_       parms_pcsetpermscaloptions__
#define parms_pcsetpermtype_              parms_pcsetpermtype__
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_pccreate_                   parms_pccreate
#define parms_pcfree_ 			  parms_pcfree
#define parms_pcgetratio_ 		  parms_pcgetratio
#define parms_pcsetbsize_ 		  parms_pcsetbsize
#define parms_pcsetfill_ 		  parms_pcsetfill
#define parms_pcsetilutype_ 		  parms_pcsetilutype
#define parms_pcsetinnereps_ 		  parms_pcsetinnereps
#define parms_pcsetinnerksize_ 		  parms_pcsetinnerksize
#define parms_pcsetinnermaxits_ 	  parms_pcsetinnermaxits
#define parms_pcsetnlevels_ 		  parms_pcsetnlevels
#define parms_pcsetparams_ 		  parms_pcsetparams
#define parms_pcsettol_ 		  parms_pcsettol
#define parms_pcsettolind_ 		  parms_pcsettolind
#define parms_pcsettype_ 		  parms_pcsettype
#define parms_pcsetup_ 			  parms_pcsetup
#define parms_pcsolve_ 			  parms_pcsolve
#define parms_pcview_            	  parms_pcview
#define parms_pcgetname_                  parms_pcgetname
#define parms_pcilugetname_               parms_pcilugetname
#define parms_pccreateabstract_           parms_pccreateabstract
#define parms_pcsetop_                    parms_pcsetop
#define parms_pcsetpermscaloptions_       parms_pcsetpermscaloptions
#define parms_pcsetpermtype_              parms_pcsetpermtype
#endif

void parms_pccreate_(parms_PC *self, parms_Mat *A, int *ierr)
{
  *ierr = parms_PCCreate(self, *A);
}

void parms_pcfree_(parms_PC *self, int *ierr)
{
  *ierr = parms_PCFree(self);
}

void parms_pcgetratio_(parms_PC *self, double *ratio, int *ierr)
{
  *ierr = parms_PCGetRatio(*self,  ratio);
}

void parms_pcsetbsize_(parms_PC *self, int *bsize, int *ierr)
{
  *ierr = parms_PCSetBsize(*self, *bsize);
}
               
void parms_pcsetfill_(parms_PC *self, int *fill, int *ierr)
{
  *ierr = parms_PCSetFill(*self, fill);
}

void parms_pcsetilutype_(parms_PC *self, PCILUTYPE *pcstype, int
			 *ierr)
{
  *ierr = parms_PCSetILUType(*self, *pcstype);
}

void parms_pcsetinnereps_(parms_PC *self, REAL *eps, int *ierr)
{
  *ierr = parms_PCSetInnerEps(*self, *eps);
}

void parms_pcsetinnerksize_(parms_PC *self, int *im, int *ierr)
{
  *ierr = parms_PCSetInnerKSize(*self, *im);
}

void parms_pcsetinnermaxits_(parms_PC *self, int *imax, int *ierr)
{
  *ierr = parms_PCSetInnerMaxits(*self, *imax);
}

void parms_pcsetnlevels_(parms_PC *self, int *nlevel, int *ierr)
{
  *ierr = parms_PCSetNlevels(*self, *nlevel);
}

void parms_pcsetparams_(parms_PC *self, int *nflags, char **params,
			int *ierr)
{
  *ierr = parms_PCSetParams(*self, *nflags, params);
}

void parms_pcsettol_(parms_PC *self, double *tol, int *ierr)
{
  *ierr = parms_PCSetTol(*self, tol);
}

void parms_pcsettolind_(parms_PC *self, REAL *tolind, int *ierr)
{
  *ierr = parms_PCSetTolInd(*self, *tolind);
}

void parms_pcsettype_(parms_PC *self, PCTYPE *pctype, int *ierr)
{
  *ierr = parms_PCSetType(*self, *pctype);
}

void parms_pcsetup_(parms_PC *self, int *ierr)
{
  *ierr = parms_PCSetup(*self);
}

void parms_pcsolve_(parms_PC *self, FLOAT *y, FLOAT *z, int *ierr)
{
  *ierr = parms_PCApply(*self, y, z);
}

void parms_pcview_(parms_PC *self, parms_Viewer *v, int *ierr)
{
  parms_PCView(*self, *v);
  *ierr = 0;
}

void parms_pcgetname_(parms_PC *self, char *name, int *size, int
		      *ierr, int len) 
{
  char *buf;
  int i;

  parms_PCGetName(*self, &buf);
  *size = strlen(buf);
  PARMS_MEMCPY(name, buf, *size);
  for (i = *size; i < len; i++) {
    name[i] = ' ';
  }

  *ierr = 0;
}

void parms_pcilugetname_(parms_PC *self, char *name, int *size, int
			 *ierr, int len)
{
  char *buf;
  int i;
  
  parms_PCILUGetName(*self, &buf);
  *size = strlen(buf); 
  PARMS_MEMCPY(name, buf, *size);
  for (i = *size; i < len; i++) {
    name[i] = ' ';
  }
  *ierr = 0;
}

void parms_pccreateabstract_(parms_PC *self, int *ierr)
{
  *ierr = parms_PCCreateAbstract(self);
}

void parms_pcsetpermscaloptions_(parms_PC *self, int *meth, int
				      *flag, int *ierr)
{
  *ierr = parms_PCSetPermScalOptions(*self, meth, *flag);
}

void parms_pcsetpermtype_(parms_PC *self, int *type, int *ierr)
{
  *ierr = parms_PCSetPermType(*self, *type);
}

void parms_pcsetop_(parms_PC *self, parms_Mat *A, int *ierr)
{
  *ierr = parms_PCSetOP(*self, *A);
}

