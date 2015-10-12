#include <ctype.h>
#include "petscksp.h"
/* REQUIRED */

#define MAX_LINE 100
#define BUFFLEN 200
#define STR2U(p, p1)				\
  if (1) {					\
    p = buf;					\
    while (' ' == *p) {				\
      p++;					\
    }						\
    p1 = malloc((strlen(p)+1)*sizeof(*p1));	\
    p2 = p1;					\
    while (' ' != *p) {				\
      *p2 = toupper(*p);			\
      ++p;					\
      ++p2;					\
    }						\
    *p2 = '\0';					\
  }
/* Define context for user-provided preconditioner */

typedef struct {
  parms_Map       map;
  parms_Mat       A;
  parms_PC        pc;
}PC_PARMS;

/* Declare routines for user-provided preconditioner */

PetscErrorCode PCCreate_PARMS(PC_PARMS **prec);
PetscErrorCode PCSetUp_PARMS(PC ctx);
PetscErrorCode PCApply_PARMS(PC dummy,Vec b,Vec x);
PetscErrorCode PCDestroy_PARMS(PC dummy);
int parms_SetOptions(parms_PC pc, char *fname);

/* 
THE FUNCTIONS BELOW ARE ONLY REQUIRED BY THE HB FILE FORMAT EXAMPLE "dd-petsc.c".
*/   
#ifndef _AUX_HEADER_INCLUDED_H
#define _AUX_HEADER_INCLUDED_H

#if defined(FORTRAN_CAPS)
#define wreadmtc_            WREADMTC
#define csrcsc_              CSRCSC
#define aplb_                APLB
#define dse_                 DSE
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define wreadmtc_            wreadmtc__
#define csrcsc_              csrcsc__
#define aplb_                aplb__
#define dse_                 dse__
#elif !defined(FORTRAN_UNDERSCORE)
#define wreadmtc_            wreadmtc
#define csrcsc_              csrcsc
#define aplb_                aplb
#define dse_                 dse
#endif

#define BUFLEN 100

/* Routine to read test matrix */
extern int read_matrix(char mname[BUFLEN]);
extern void wreadmtc_(int *nmax, int *nzmax, int *job, char *fname,int
		      *len, double *a, int *ja, int *ia, double *rhs,
		      int *nrhs, char *guesol, int *nrow, int *ncol,
		      int *nnz, char *title, char *key, char *type,
		      int *ierr); 
extern void csrcsc_(int *n, int *job, int *ipos, double *a, int *ja,
		    int *ia, double *ao, int *jao, int *iao);  
extern void aplb_(int *nrow, int *ncol, int *job, double *a, int *ja,
		  int *ia, double *b, int *jb, int *ib, double *c, int
		  *jc, int *ic, int *nnzmax, int *iw, int *ierr); 



/*---------- complex routines ------------*/
#if defined(DBL_CMPLX)
#if defined(FORTRAN_CAPS)
#define wreadmtc_            WREADMTC
#define readmtc_             READMTC
#define csrcsc_              CSRCSC
#define aplb_                APLB
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define wreadmtc_            wreadmtc__
#define readmtc_             readmtc__
#define csrcsc_              csrcsc__
#define aplb_                aplb__
#elif !defined(FORTRAN_UNDERSCORE)
#define wreadmtc_            wreadmtc
#define readmtc_             readmtc
#define csrcsc_              csrcsc
#define aplb_                aplb
#endif
extern void zwreadmtc_(int *nmax, int *nzmax, int *job, char *fname,int
		      *len, complex double *a, int *ja, int *ia, complex double *rhs,
		      int *nrhs, char *guesol, int *nrow, int *ncol,
		      int *nnz, char *title, char *key, char *type,
		      int *ierr); 
extern void zreadmtc_(int *nmax, int *nzmax, int *job, char *fname, complex double 
		     *a, int *ja, int *ia, complex double *rhs, int *nrhs, char 
		     *guesol, int *nrow, int *ncol,int *nnz, char *title, 
		     char *key, char *type, int *ierr);
extern void zcsrcsc_(int *n, int *job, int *ipos, complex double *a, int *ja,
		    int *ia, complex double *ao, int *jao, int *iao);  
extern void zaplb_(int *nrow, int *ncol, int *job, complex double *a, int *ja,
		  int *ia, complex double *b, int *jb, int *ib, complex double *c, int
		  *jc, int *ic, int *nnzmax, int *iw, int *ierr); 
#endif

#endif 

