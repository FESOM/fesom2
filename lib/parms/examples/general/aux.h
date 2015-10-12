#ifndef _AUX_HEADER_INCLUDED_H
#define _AUX_HEADER_INCLUDED_H
#define  MAX_LINE  100
#define MAX_MAT 100

#include "parms.h"

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
  } else

typedef struct fprm_ {
  PCTYPE     pctype;
  PCILUTYPE  pcilutype;
  double     pgfpar[2];
  double     tolind;
  double     droptol[7];
  int        ipar[18];
  int        lfil[7];
} *fprm;

#if defined(FORTRAN_CAPS)
#define fread_param_         FREAD_PARAM
#define fset_pc_params_      FSET_PC_PARAMS
#define fset_solver_params_  FSET_SOLVER_PARAMS
#define fprm_free_           FPRM_FREE
#define wreadmtc_            WREADMTC
#define readmtc_             READMTC
#define csrcsc_              CSRCSC
#define aplb_                APLB
#define dse_                 DSE
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define fread_param_         fread_param__
#define fset_pc_params_      fset_pc_params__
#define fset_solver_params_  fset_solver_params__
#define fprm_free_           fprm_free__
#define wreadmtc_            wreadmtc__
#define readmtc_             readmtc__
#define csrcsc_              csrcsc__
#define aplb_                aplb__
#define dse_                 dse__
#elif !defined(FORTRAN_UNDERSCORE)
#define fread_param_         fread_param
#define fset_pc_params_      fset_pc_params
#define fset_solver_params_  fset_solver_params
#define fprm_free_           fprm_free
#define wreadmtc_            wreadmtc
#define readmtc_             readmtc
#define csrcsc_              csrcsc
#define aplb_                aplb
#define dse_                 dse
#endif

extern int  read_param(char *fname, char mname[MAX_MAT][MAX_LINE], fprm prm);
extern void set_pc_params(parms_PC pc, fprm prm);
extern void set_solver_params(parms_Solver solver, fprm prm);
extern void fread_param_(char *fname, fprm *prm, char *matrix, int *matlen, int len);
extern void fset_pc_params_(parms_PC *pc, fprm *prm);
extern void fset_solver_params_(parms_Solver *solver, fprm *prm);
extern void fprm_free_(fprm *prm);
extern void wreadmtc_(int *nmax, int *nzmax, int *job, char *fname,int
		      *len, double *a, int *ja, int *ia, double *rhs,
		      int *nrhs, char *guesol, int *nrow, int *ncol,
		      int *nnz, char *title, char *key, char *type,
		      int *ierr); 
extern void readmtc_(int *nmax, int *nzmax, int *job, char *fname, double 
		     *a, int *ja, int *ia, double *rhs, int *nrhs, char 
		     *guesol, int *nrow, int *ncol,int *nnz, char *title, 
		     char *key, char *type, int *ierr);
extern void csrcsc_(int *n, int *job, int *ipos, double *a, int *ja,
		    int *ia, double *ao, int *jao, int *iao);  
extern void aplb_(int *nrow, int *ncol, int *job, double *a, int *ja,
		  int *ia, double *b, int *jb, int *ib, double *c, int
		  *jc, int *ic, int *nnzmax, int *iw, int *ierr); 
extern void dse_(int *n, int *ja, int *ia, int *ndom, int *riord, int
		 *dom, int *idom, int *mask, int *jwk, int *link);

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
