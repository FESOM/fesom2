#ifndef PARMS_AUX_READ_H
#define PARMS_AUX_READ_H

#include <string.h>
#include <stdlib.h>
#include "parms.h"

#define MAXLINE 200
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
  int        mprocx;
  int        mprocy;
  int        xnmesh;
  int        ynmesh;
} *fprm;

#if defined(FORTRAN_CAPS)
#define fread_param_         FREAD_PARAM
#define fset_pc_params_      FSET_PC_PARAMS
#define fset_solver_params_  FSET_SOLVER_PARAMS
#define fprm_free_           FPRM_FREE
#define part1_               PART1
#define part2_               PART2
#define gen5loc_             GEN5LOC
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define fread_param_         fread_param__
#define fset_pc_params_      fset_pc_params__
#define fset_solver_params_  fset_solver_params__
#define fprm_free_           fprm_free__
#define part1_               part1__
#define part2_               part2__
#define gen5loc_             gen5loc__
#elif !defined(FORTRAN_UNDERSCORE)
#define fread_param_         fread_param
#define fset_pc_params_      fset_pc_params
#define fset_solver_params_  fset_solver_params
#define fprm_free_           fprm_free
#define part1_               part1
#define part2_               part2
#define gen5loc_             gen5loc
#endif

extern int  read_param(char *fname, fprm prm);
extern void set_pc_params(parms_PC pc, fprm prm);
extern void set_solver_params(parms_Solver solver, fprm prm);
extern void fset_pc_params_(parms_PC *pc, fprm *prm);
extern void fset_solver_params_(parms_Solver *solver, fprm *prm);
extern void fread_param_(char *fname, fprm *prm, int *mprocx, int
			 *mprocy, int *xnmesh, int *ynmesh, int len);
extern void fprm_free_(fprm *prm);
extern void part1_(int *nx,int *ny,int *mpx,int *mpy,int *ovp,int
		   *lst,int *lstptr,int *iout); 
extern void part2_(int *nx,int *ny,int *nz,int *mpx,int *mpy,int *mpz,
		   int *ovp,int *lst,int *lstptr,int *iout); 
extern void gen5loc_(int *nx,int *ny,int *nz,int *nloc,int
		     *riord,double *a,int *ja,int *ia, double
		     *stencil); 
#if defined(DBL_CMPLX)		     
extern void zgen5loc_(int *nx,int *ny,int *nz,int *nloc,int
		     *riord,complex double *a,int *ja,int *ia, complex double
		     *stencil);
#endif		     
#endif 
