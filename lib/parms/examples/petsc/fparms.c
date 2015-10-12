/* This file contains pARMS interface routines to petsc for fortran calls. Note that this is 
 * only compatible with the latest version of petsc (version 3.1). For older 
 * versions (version 3.0) use fparms.c_petsc_3.0.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "petsc.h"
#include "petscksp.h"
#include "petscsys.h"
#include "parms.h"
#include "protos.h"
#ifdef  _SUN 
#define f_pcsetup_parms_ f_pcsetup_parms_
#define f_pcapply_parms_ f_pcapply_parms_
#define f_pcdestroy_parms_ f_pcdestroy_parms_
#else 
#ifdef _LINUX
#define f_pcsetup_parms_ f_pcsetup_parms__
#define f_pcapply_parms_ f_pcapply_parms__
#define f_pcdestroy_parms_ f_pcdestroy_parms__
#else 
#ifdef _IBM
#define f_pcsetup_parms_ f_pcsetup_parms
#define f_pcapply_parms_ f_pcapply_parms
#define f_pcdestroy_parms_ f_pcdestroy_parms
#else 
#ifdef _SGI 
#define f_pcsetup_parms_ f_pcsetup_parms_
#define f_pcapply_parms_ f_pcapply_parms_
#define f_pcdestroy_parms_ f_pcdestroy_parms_
#endif
#endif
#endif
#endif

/*---------------------------------------------------------------------
| interface  for FORTRAN calls to pARMS routines. 
+---------------------------------------------------------------------*/

PetscErrorCode f_pcsetup_parms_(PC *ctx, PetscErrorCode *ierr)
{
/* function protos */
  PetscErrorCode PCSetUp_PARMS(PC );
  int parms_SetOptions(parms_PC , char *);
/*end protos */

  *ierr = PCSetUp_PARMS(*ctx);CHKERRQ(*ierr);

  return(*ierr);

}

PetscErrorCode f_pcapply_parms_(PC *dummy, Vec *b,Vec *x, PetscErrorCode *ierr)
{
/*function protos */
  PetscErrorCode PCApply_PARMS(PC,Vec ,Vec );
/*end protos*/

  *ierr = PCApply_PARMS(*dummy,*b,*x);CHKERRQ(*ierr);

  return(*ierr);
}

PetscErrorCode f_pcdestroy_parms_(PC *dummy, PetscErrorCode *ierr)
{

/* function protos */
  PetscErrorCode PCDestroy_PARMS(PC);
/* end protos */

  *ierr = PCDestroy_PARMS(*dummy);CHKERRQ(*ierr);

  return(*ierr);
}


