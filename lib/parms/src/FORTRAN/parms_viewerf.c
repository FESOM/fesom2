#include "parms_mem.h"
#include "parms_viewer.h"

#if defined(FORTRAN_CAPS)
#define parms_viewercreate_      PARMS_VIEWERCREATE
#define parms_viewerfree_	 PARMS_VIEWERFREE
#define parms_viewergetfp_	 PARMS_VIEWERGETFP
#define parms_viewergetfname_	 PARMS_VIEWERGETFNAME
#define parms_viewerstorefp_     PARMS_VIEWERSTOREFP
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_viewercreate_      parms_viewercreate__
#define parms_viewerfree_	 parms_viewerfree__
#define parms_viewergetfp_	 parms_viewergetfp__
#define parms_viewergetfname_	 parms_viewergetfname__
#define parms_viewerstorefp_     parms_viewerstorefp__
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_viewercreate_      parms_viewercreate
#define parms_viewerfree_	 parms_viewerfree
#define parms_viewergetfp_	 parms_viewergetfp
#define parms_viewergetfname_	 parms_viewergetfname
#define parms_viewerstorefp_     parms_viewerstorefp
#endif

void parms_viewercreate_(parms_Viewer *self, char *fname, int *ierr,
			 int len)
{
  char *buf;

  PARMS_NEWARRAY(buf, len+1);
  PARMS_MEMCPY(buf, fname, len);
  buf[len] = '\0';
  *ierr = parms_ViewerCreate(self, buf);
  PARMS_FREE(buf);
}

void parms_viewerfree_(parms_Viewer *self, int *ierr)
{
  *ierr = parms_ViewerFree(self);
}

void parms_viewergetfp_(parms_Viewer *self, FILE **fp, int *ierr)
{
  *ierr = parms_ViewerGetFP(*self, fp);
}

void parms_viewergetfname_(parms_Viewer *self, char **fname, int
			   *ierr)
{
  *ierr = parms_ViewerGetFname(*self, fname);
}

void parms_viewerstorefp_(parms_Viewer *self, FILE *fp, int *ierr)
{
  *ierr = parms_ViewerStoreFP(*self,  fp);
}
