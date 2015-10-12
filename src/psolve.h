#ifndef __PSOLVE__
#define __PSOLVE__

#include "parms.h"
#include "parms_sys.h"

typedef struct psolver_{

    int reuse;
    
    int *rows;
    int *rptr;
    int *cols;
    double *vals;
    double *scale;
    parms_Solver ksp;
    parms_Map map;

} *psolver;

#ifdef FORTRAN_NOUNDERSCORE
#define psolver_init_ psolver_init
#define psolve_ psolve
#define psolver_final_ psolver_final
#endif

int set_pc_params(parms_PC pc, PCTYPE pctype, PCILUTYPE pcilutype, 
		  int ilulevel, int fillin, double droptol);

int set_solver_params(parms_Solver solver, SOLVERTYPE solvertype, 
		      int maxits, int restart, double soltol);

extern void psolver_init_(int *id, SOLVERTYPE *stype, PCTYPE *pctype, PCILUTYPE *pcilutype,
		          int *ilulevel, int *fillin, double *droptol, int *maxits, int *restart, double *soltol, 
		          int *part, int *rptr, int *cols, double *vals, int *reuse);
extern void psolve_(int *id, double *rhs, double *vals, double *sol, int *new);
extern void psolver_final_();

#endif
