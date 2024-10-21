#ifndef __PSOLVE__
#define __PSOLVE__

#include "parms.h"
#include "parms_sys.h"

typedef struct psolver{

    int reuse;
    
    int *rows;
    int *rptr;
    int *cols;
    double *vals;
    double *scale;
    parms_Solver ksp;
    parms_Map map;

} *psolver;

int set_pc_params(parms_PC pc, PCTYPE pctype, PCILUTYPE pcilutype, 
		  int ilulevel, int fillin, double droptol);

int set_solver_params(parms_Solver solver, SOLVERTYPE solvertype, 
		      int maxits, int restart, double soltol);

#endif
