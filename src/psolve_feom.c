#ifdef PARMS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "psolve.h"

#define NSOL 1 
//#define NSOL 10

psolver solvers[NSOL];
int solv_id[12] = {0};
int nsolver = 0;

void psolver_init_(int *id, SOLVERTYPE *stype, PCTYPE *pctype, PCILUTYPE *pcilutype,
		   int *ilulevel, int *fillin, double *droptol, int *maxits, int *restart, double *soltol, 
		   int *part, int *rptr, int *cols, double *vals, int *reuse)
{
  
  parms_Viewer v;
  int i, j, k, nloc, pid, nproc;
  int *ncnts, *idxn, *rp=NULL, *r=NULL, *c=NULL, nmb;
  double  tmp, *scale, *values=NULL;
 
 parms_Map map;
  parms_Mat A;
  parms_PC pc;
  parms_Solver ksp;

  psolver solver;
  
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  
  solver = malloc(sizeof(*solver));
  nmb = part[nproc]-part[0];

  nloc = 0;
  if(nproc > 1) {
      nloc = part[pid+1]-part[pid];
      parms_MapCreateFromPetsc(&map, nloc, nmb, MPI_COMM_WORLD);        
      idxn = malloc(nloc*sizeof(int));    
      parms_MapGetGlobalIndices(map, idxn);
  }
  else {
      parms_MapCreateFromLocal(&map,nmb,0);
      nloc = nmb;
      idxn = malloc(nloc*sizeof(int));    
      for(i = 0; i < nloc; i++)
	  idxn[i] = i;
  }    
  solver->reuse = *reuse;
  
  scale  = malloc(nloc*sizeof(double));
  values = malloc(rptr[nloc]*sizeof(double));
  for(i = 0; i < nloc; i++){
      tmp = 0.;
      for(j = rptr[i]; j < rptr[i+1]; j++)
	tmp += fabs(vals[j]);
      scale[i] = 1./tmp;
      for(j = rptr[i]; j < rptr[i+1]; j++)
	values[j] = vals[j]*scale[i];
  }
  solver->scale = scale;
  
  // create Mat
  parms_MatCreate(&A,map);
  parms_MatSetValues(A, nloc, idxn, rptr, cols, values, INSERT);
  parms_MatSetup(A);
  
  // create PC/Solver & set parameters
  parms_PCCreate(&pc,A);
  set_pc_params(pc, *pctype, *pcilutype, *ilulevel, *fillin, *droptol);
  parms_PCSetup(pc);
  parms_SolverCreate(&ksp,A,pc);
  set_solver_params(ksp, *stype, *maxits, *restart, *soltol); 

  solver->ksp = ksp;
  solver->map = map; 
 
  if(solver->reuse){
      r = malloc((nloc)*sizeof(int));
      for(i = 0; i < nloc; i++)
        r[i] = idxn[i];    
      solver->rows = r;

  rp = malloc((nloc+1)*sizeof(int));
  for(i = 0; i < nloc+1; i++)
      rp[i] = rptr[i];
      solver->rptr = rp;
  
  c = malloc(rptr[nloc]*sizeof(int));
  for(i = 0; i < rptr[nloc]; i++)
    c[i] = cols[i];
  solver->cols = c;
  solver->vals = values;
  }
  else{
     if(!r)
        free(r);
     if(!rp)
        free(rp);
     if(!c)
        free(c);
     if(!values)
        free(values);
  }
      

  solv_id[*id] = nsolver;
  solvers[nsolver++] = solver;
   
}

void psolver_final_()
{

  int pid, i;
  parms_Solver ksp;
  parms_PC pc;
  parms_Mat A;
  parms_Map map;
  psolver solver;

  for(i = 0; i < nsolver; i++){
      solver = solvers[i];
      
      ksp = solver->ksp;
      parms_SolverGetPC(ksp, &pc);
      parms_SolverGetMat(ksp,&A);
      map = solver->map;
      
      parms_SolverFree(&ksp);
      parms_PCFree(&pc);
      parms_MatFree(&A);
      parms_MapFree(&map);
 
      if(solver->reuse){
	  free(solver->rptr);
	  free(solver->rows);
	  free(solver->cols);
	  free(solver->vals);
      }
      free(solver->scale);

      free(solver);
  }
}

void psolve_(int *id, double *rhs, double *vals, double *sol, int *new)
{
  
  parms_Viewer v;
  psolver solver;
  parms_Map map;
  parms_Mat A;
  parms_PC pc;
  parms_Solver ksp;

  double resnorm;
  double *x,*y, *scale, *values, tmp;
  double t0,t1,t2,toh;
  int *rptr, *cols;
  int its,err,i,j,k, cnt, nloc, pid, sid;
  
  sid = solv_id[*id];

  solver = solvers[sid];
  ksp = solver->ksp;
  map = solver->map;
  nloc = parms_MapGetLocalSize(map);
  
  scale = solver->scale;
 
  if(*new) {
    if(solver->reuse){
      parms_SolverGetPC(ksp, &pc);
      parms_SolverGetMat(ksp, &A);
    
      rptr = solver->rptr;
      cols = solver->cols;
    
      values = solver->vals;
      for(i = 0; i < nloc; i++){
	tmp = 0.;
	for(j = rptr[i]; j < rptr[i+1]; j++)
	  tmp += fabs(vals[j]);
	scale[i] = 1./tmp;
	for(j = rptr[i]; j < rptr[i+1]; j++)
	  values[j] = vals[j]*scale[i];
      }
    
      // create Mat & set values
      parms_MatReset(A,SAME_NONZERO_STRUCTURE); 
      parms_MatSetValues(A, nloc, solver->rows, rptr, cols, values, INSERT);
      parms_MatSetup(A); 
    
      parms_PCSetup(pc);
    }
    else
      printf("ERROR: matrix data is static\n");
  }

  x = sol;
  y = malloc(nloc*sizeof(double));
  for(i = 0; i < nloc; i++)
    y[i] = rhs[i]*scale[i];  	  
  
  // solve system of equations  
  parms_SolverApply(ksp,y,x);
  
/*  
  // get trueresidual and number of iterations
  parms_SolverGetResidualNorm2(ksp,y,x,&resnorm);
  its = parms_SolverGetIts(ksp);
  printf("%e %d\n", resnorm, its);
*/  

  free(y);
}

int set_pc_params(parms_PC pc, PCTYPE pctype, PCILUTYPE pcilutype, 
		  int ilulevel, int fillin, double droptol){

  int i, lfil[7];
  double dtol[7];
  
  for(i = 0; i < 7; i++){
    lfil[i] = fillin;
    dtol[i] = droptol;
  }
    
  parms_PCSetType(pc,    pctype);
  parms_PCSetILUType(pc, pcilutype);
  parms_PCSetNlevels(pc, ilulevel);
  parms_PCSetFill(pc,    lfil);
  parms_PCSetTol(pc,     dtol);
  
  return 0;
}


int set_solver_params(parms_Solver solver, SOLVERTYPE solvertype, 
		      int maxits, int restart, double soltol){
    char buf[100];
    
    parms_SolverSetType(solver, solvertype);
  
    sprintf(buf, "%d", maxits);
    parms_SolverSetParam(solver, MAXITS, buf);

    sprintf(buf, "%d", restart);
    parms_SolverSetParam(solver, KSIZE,  buf);
    
    sprintf(buf, "%g", soltol);
    parms_SolverSetParam(solver, DTOL,  buf);
    
    return 0;
}


#endif
