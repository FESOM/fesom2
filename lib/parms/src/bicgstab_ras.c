/*--------------------------------------------------------------------
  bicgstabras_create    : create the BICGS_RAS solver. 
  bicgstabras_free      : free the memory for the BICGS_RAS solver. 
  bicgstabras_view      : dump the BICGS_RAS solver.
  parms_bicgstabras     : the BICGS_RAS solve function.

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

#include "./include/parms_comm_impl.h"
#include "./include/parms_pc_impl.h"
#include "./include/parms_opt_impl.h"


typedef struct ras_data {
  parms_Operator op;
  parms_Comm handler;
  parms_Map  is;
  BOOL issetup;
  FLOAT *rbuf;
  int nloc;
  int n;
  int nodv;
  int nsend;
} *ras_data;


#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

typedef struct bicgs_data {
  int neigs;
} *bicgstabras_data;

#define TINY 1.0e-20
int parms_bicgstabras(parms_Solver self, FLOAT *y, FLOAT *x)
{

    int i, its;
    int maxits, nloc, n_ext, size, one = 1;
    FLOAT alpha, t1, tmp[2];
    FLOAT *r0, *v, *t;
    FLOAT *r_ext, *st_ext, *p_ext, *pt_ext;
    REAL  tol, rho, rho_alt, omega, sigma, beta;
    MPI_Request req;
    int nodv, nsend;

    parms_Mat A;
    parms_PC  pc;
    parms_Map is;
    parms_Viewer viewer;
    parms_Comm pc_handler;
    ras_data pc_data;
    
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
    
    // For inlined RAS preconditioner
    pc_data = (ras_data)pc->data;
    n_ext = pc_data->n;  // local size with halo for Schwarz
    pc_handler = pc_data->handler;
    nsend = pc_data->nsend;
    nodv = pc_data->nodv;

/*----- allocate memory for working local arrays -------*/
    size = nloc;
    /*    PARMS_NEWARRAY(r0, size);  */
    PARMS_NEWARRAY(r0,  size);   
    PARMS_NEWARRAY(t, size);
    PARMS_NEWARRAY(v, size);

    PARMS_NEWARRAY(r_ext,  n_ext); 
    PARMS_NEWARRAY(p_ext,  n_ext);
    PARMS_NEWARRAY(pt_ext, n_ext);
    PARMS_NEWARRAY(st_ext, n_ext);


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
    parms_MatVec(A, x, r_ext);
    parms_VecAYPX(r_ext, r0, -1.0, is);    
    
    /* choose r0 (arbitrary), e.g., PARMS_MEMCPY(r0,r,nloc); or: */    
    /* PARMS_MEMCPY(r0, y, nloc); */
    /* As r0 and y are never changed, we can simply set */
    /* r0=y; */
    /* Do this already above, use r0 as permuted y! */

    tmp[0] = 0.;
    tmp[1] = 0.;
    for (i = 0; i < nloc; i++) 
      {
	tmp[0] +=  r_ext[i] * r_ext[i]; // omega
	tmp[1] += r0[i] * r_ext[i]; // rho
      }
     if(is->isserial == false)
       MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req); 
     

    PARMS_MEMCPY(p_ext, r_ext, nloc);
    alpha = 0.;
    sigma = rho_alt = 1.0;

    if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);
    
    omega = tmp[0];
    rho   = tmp[1];
   
 
    if (omega > tol){
      for (its = 0; its < maxits; its++){
                 
	// parms_PCApply(pc, p, pt);
	if (nsend) parms_CommDataBegin(pc_handler, p_ext, 0);
	if (nodv) parms_CommDataEnd(pc_handler);
	if (pc_data->n - pc_data->nloc) 
	  PARMS_MEMCPY(&p_ext[nloc], pc_data->rbuf, pc_data->n - pc_data->nloc);
	
	/* solve the extended linear system */
	parms_OperatorApply(pc_data->op, p_ext, pt_ext);
      
      /* No halo exchange of the solution, that's the "Restricted" in RAS */
      /*  === End Preconditioner === */

	parms_MatVec(A, pt_ext, v);
	parms_VecDOT(r0,v,&alpha,is);

	alpha = rho/alpha;

	for (i = 0; i < nloc; i++) r_ext[i] -= alpha * v[i]; 
      
	// parms_PCApply(pc, r, st);
	if (nsend) parms_CommDataBegin(pc_handler, r_ext, 0);
	if (nodv) parms_CommDataEnd(pc_handler);
	if (pc_data->n - pc_data->nloc) 
	  PARMS_MEMCPY(&r_ext[nloc], pc_data->rbuf, pc_data->n - pc_data->nloc);
	
	/* solve the extended linear system */
	parms_OperatorApply(pc_data->op, r_ext, st_ext);
      
      /* No halo exchange of the solution, that's the "Restricted" in RAS */
      /*  === End Preconditioner === */
	parms_MatVec(A, st_ext, t);

	tmp[0] = 0.;
	tmp[1] = 0.;
	for (i = 0; i < nloc; i++) 
	  {
	    tmp[0] += t[i] * r_ext[i]; 
	    tmp[1] += t[i] * t[i]; 
	  }
            
	
	if(is->isserial == false) 
	  MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req);
    	
	//for (i = 0; i < nloc; i++) x[i] += alpha * pt_ext[i];
    
	if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);

	sigma = tmp[0]/tmp[1];
      
      
	for (i = 0; i < nloc; i++) r_ext[i] -= sigma * t[i];
     
	tmp[0] = 0.;
	tmp[1] = 0.;
	for (i = 0; i < nloc; i++) 
	  {
	    tmp[0] +=  r_ext[i] * r_ext[i]; // omega
	    tmp[1] += r0[i] * r_ext[i]; // rho
	  }
	
	if(is->isserial == false) 
	  MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req);

	for (i = 0; i < nloc; i++) x[i] += alpha * pt_ext[i] + sigma * st_ext[i];

	
	 if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);

	omega = tmp[0];

	if (omega < tol) break;

	rho   = tmp[1];

	beta = (rho * alpha) / (rho_alt * sigma);
	rho_alt = rho;
      
	for (i = 0; i < nloc; i++) 
	  {
	    p_ext[i] = r_ext[i] + beta*(p_ext[i] - sigma* v[i]);
	  }
      }
    }
    

    its++;

    if(its == maxits && omega > tol)
      printf("ERROR: no convergence\n");
    
    /* reset isvecperm and do inverse permutation*/
    if (is->isperm) { 
        for (i = 0; i < nloc; i++) {
	  t[is->iperm[i]] = x[i];
	}
	memcpy(x, t, nloc*sizeof(FLOAT));
	is->isvecperm = false; 
    }
    free(r0); 
    free(r_ext);
    free(p_ext);
    free(st_ext);
    free(t);
    free(v);
    free(pt_ext);
    

    self->its = its;
    return 0;
}

static int bicgstabras_getresidual(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *res)
{

    parms_Mat A;

    A  = self->A;

    parms_MatMVPY(A, -1.0, x, 1.0, y, res);

    return 0;
}
     
static int bicgstabras_getresidualnorm2(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm)
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


static int bicgstabras_free(parms_Solver *self)
{
  /*dummy*/
  return 0;
}

static int bicgstabras_view(parms_Solver self, parms_Viewer v)
{
  FILE *fp;
  char *name, *iluname;
  int dim;
  parms_ViewerGetFP(v, &fp);

  fprintf(fp,"\n=============================\n");
  fprintf(fp,"	Solver Parameters	\n");
  fprintf(fp,"=============================\n");
  
  fprintf(fp, "Solver type = BiConjugate Gradient Stabilized Method (BICGS_RAS) \n");
  
  fprintf(fp, "maxits = %d \n", self->maxits);
  fprintf(fp, "Relative tolerance = %-8.2e \n", self->tol);

  parms_PCGetName(self->pc, &name);
  parms_PCILUGetName(self->pc, &iluname);
  fprintf(fp, "Global Preconditioner: %s\n", name);
  fprintf(fp, "Local Preconditioner: %s\n", iluname);

  parms_ViewerGetFP(v, &fp);

  return 0;
}

static int bicgstabras_setksize(parms_Solver self, int restart)
{
  /*dummy*/
  return 0;
}

static int bicgstabras_setneig(parms_Solver self, int neigs)
{
  /*dummy*/
  return 0;
}

static struct parms_Solver_ops parms_bicgstabras_sol = {
  parms_bicgstabras,
  bicgstabras_getresidual,
  bicgstabras_getresidualnorm2,
  bicgstabras_setksize,
  bicgstabras_setneig,
  bicgstabras_free,
  bicgstabras_view
};


/** Create the BICGS_RAS solver. 
 *
 *  \param self A parms_Solver object.
 *  \return 0 on success.
 */
int bicgstabras_create(parms_Solver self)
{
  PARMS_MEMCPY(self->ops, &parms_bicgstabras_sol, 1);		
  return 0;
}
