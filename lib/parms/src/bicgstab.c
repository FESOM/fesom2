/*--------------------------------------------------------------------
  bicgstab_create    : create the BICGS solver. 
  bicgstab_free      : free the memory for the BICGS solver. 
  bicgstab_view      : dump the BICGS solver.
  parms_bicgstab     : the BICGS solve function.

  $Id: fgmres.c,v 1.4 2006-12-15 07:02:07 zzli Exp $
  ------------------------------------------------------------------*/

#if defined(__ICC)
#include <mathimf.h>
#else
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#endif 
#include "parms_vec.h"
#include "parms_mat.h"
#include "parms_pc.h"
#include "./include/parms_solver_impl.h"

#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

typedef struct bicgs_data {
  int neigs;
} *bicgstab_data;

#define TINY 1.0e-20
int parms_bicgstab(parms_Solver self, FLOAT *y, FLOAT *x)
{

    int i, its;
    int maxits, nloc, size, one = 1;
    FLOAT alpha, t1, tmp[2];
    FLOAT *r0, *r, *p, *st, *v, *t, *pt;
    REAL  tol, rho, rho_alt, omega, sigma, beta;
    MPI_Request req;

    parms_Mat A;
    parms_PC  pc;
    parms_Map is;
    parms_Viewer viewer;
    
    int rank;
    MPI_Comm comm;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
    A  = self->A;    
    pc = self->pc;
    is = A->is;
   
    maxits  = self->maxits;
    tol     = self->tol*self->tol;
    nloc    = parms_MapGetLocalSize(is);
    comm    = is->comm;
    
/*----- allocate memory for working local arrays -------*/
    size = nloc;
    /*    PARMS_NEWARRAY(r0, size);  */
    PARMS_NEWARRAY(r,  size); 
    PARMS_NEWARRAY(r0,  size);   
    PARMS_NEWARRAY(st, size);
    PARMS_NEWARRAY(t, size);
    PARMS_NEWARRAY(v, size);
    PARMS_NEWARRAY(p,  size);
    PARMS_NEWARRAY(pt, size);

/* --- first permute x and y for pARMS matrix structure ------*/
    if (is->isperm) { 
      for (i = 0; i < nloc; i++) {
	t[is->perm[i]] = x[i];
	r0[is->perm[i]] = y[i];
      }
      memcpy(x, t, nloc*sizeof(FLOAT));
      is->isvecperm = true;
    }
   
   
    
    /* compute residual vector r = y - Ax */
    /* r0 = y  (permuted) */
    parms_MatVec(A, x, r);
    parms_VecAYPX(r, r0, -1.0, is);    
    
    /* choose r0 (arbitrary), e.g., PARMS_MEMCPY(r0,r,nloc); or: */    
    /* PARMS_MEMCPY(r0, y, nloc); */
    /* As r0 and y are never changed, we can simply set */
    /* r0=y; */
    /* Do this already above, use r0 as permuted y! */

    tmp[0] = 0.;
    tmp[1] = 0.;
    for (i = 0; i < nloc; i++) 
      {
	tmp[0] +=  r[i] * r[i]; // omega
	tmp[1] += r0[i] * r[i]; // rho
      }
     if(is->isserial == false)
       MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req); 
     

    PARMS_MEMCPY(p, r, nloc);
    alpha = 0.;
    sigma = rho_alt = 1.0;

    if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);
    
    omega = tmp[0];
    rho   = tmp[1];
   
 
    if (omega > tol){
      for (its = 0; its < maxits; its++){                 

	parms_PCApply(pc, p, pt);
	parms_MatVec(A, pt, v);
	parms_VecDOT(r0,v,&alpha,is);

	alpha = rho/alpha;

	for (i = 0; i < nloc; i++) r[i] -= alpha * v[i]; 
      
	parms_PCApply(pc, r, st);
	parms_MatVec(A, st, t);

	tmp[0] = 0.;
	tmp[1] = 0.;
	for (i = 0; i < nloc; i++) 
	  {
	    tmp[0] += t[i] * r[i]; 
	    tmp[1] += t[i] * t[i]; 
	  }
            
	
	if(is->isserial == false) 
	  MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req);
    	
	for (i = 0; i < nloc; i++) x[i] += alpha * pt[i];
    
	if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);

	sigma = tmp[0]/tmp[1];
      
	if(ABS_VALUE(sigma) < TINY){
	  if(rank == 0)
	    printf("ERROR: sigma = 0\n");
	  break;
	}
      
	for (i = 0; i < nloc; i++) r[i] -= sigma * t[i];
     
	tmp[0] = 0.;
	tmp[1] = 0.;
	for (i = 0; i < nloc; i++) 
	  {
	    tmp[0] +=  r[i] * r[i]; // omega
	    tmp[1] += r0[i] * r[i]; // rho
	  }
	
	if(is->isserial == false) 
	  MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req);

	for (i = 0; i < nloc; i++) x[i] += sigma * st[i];

	
	 if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);

	omega = tmp[0];

	if (omega < tol) break;

	rho   = tmp[1];

	beta = (rho * alpha) / (rho_alt * sigma);
	rho_alt = rho;
      
	for (i = 0; i < nloc; i++) 
	  {
	    p[i] = r[i] + beta*(p[i] - sigma* v[i]);
	  }
      }
    }
    

    its++;

    if(its == maxits && omega > tol)
      printf("ERROR: no convergence\n");
    
    /* reset isvecperm and do inverse permutation*/
    if (is->isperm) { 
        for (i = 0; i < nloc; i++) {
	  r[is->iperm[i]] = x[i];
	}
	memcpy(x, r, nloc*sizeof(FLOAT));
	is->isvecperm = false; 
    }
    free(r0); 
    free(r);
    free(p);
    free(st);
    free(t);
    free(v);
    free(pt);
    

    self->its = its;
    return 0;
}

static int bicgstab_getresidual(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *res)
{

    parms_Mat A;

    A  = self->A;

    parms_MatMVPY(A, -1.0, x, 1.0, y, res);

    return 0;
}
     
static int bicgstab_getresidualnorm2(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm)
{

    int nloc;
    FLOAT *res;
    parms_Mat A;
    parms_Map is;

    A  = self->A;
    is = A->is;
    nloc = parms_MapGetLocalSize(is);
    PARMS_NEWARRAY(res,  nloc);    

    parms_MatMVPY(A, -1.0, x, 1.0, y, res);    

    parms_VecGetNorm2(res, rnorm, is);

    PARMS_FREE(res);

    return 0;
}    


static int bicgstab_free(parms_Solver *self)
{
  /*dummy*/
  return 0;
}

static int bicgstab_view(parms_Solver self, parms_Viewer v)
{
  FILE *fp;
  char *name, *iluname;
  int dim;
  parms_ViewerGetFP(v, &fp);

  fprintf(fp,"\n=============================\n");
  fprintf(fp,"	Solver Parameters	\n");
  fprintf(fp,"=============================\n");
  
  fprintf(fp, "Solver type = BiConjugate Gradient Stabilized Method (BICGS) \n");
  
  fprintf(fp, "maxits = %d \n", self->maxits);
  fprintf(fp, "Relative tolerance = %-8.2e \n", self->tol);

  parms_PCGetName(self->pc, &name);
  parms_PCILUGetName(self->pc, &iluname);
  fprintf(fp, "Global Preconditioner: %s\n", name);
  fprintf(fp, "Local Preconditioner: %s\n", iluname);

  parms_ViewerGetFP(v, &fp);

  return 0;
}

static int bicgstab_setksize(parms_Solver self, int restart)
{
  /*dummy*/
  return 0;
}

static int bicgstab_setneig(parms_Solver self, int neigs)
{
  /*dummy*/
  return 0;
}

static struct parms_Solver_ops parms_bicgstab_sol = {
  parms_bicgstab,
  bicgstab_getresidual,
  bicgstab_getresidualnorm2,
  bicgstab_setksize,
  bicgstab_setneig,
  bicgstab_free,
  bicgstab_view
};


/** Create the BICGS solver. 
 *
 *  \param self A parms_Solver object.
 *  \return 0 on success.
 */
int bicgstab_create(parms_Solver self)
{
  PARMS_MEMCPY(self->ops, &parms_bicgstab_sol, 1);		
  return 0;
}
