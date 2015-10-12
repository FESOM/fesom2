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
    FLOAT *r0, *r, *p, *s, *st, *v, *t, *pt;
    REAL  tol, rho, rho_alt, omega, sigma, beta;

    parms_Mat A;
    parms_PC  pc;
    parms_Map is;
    parms_Viewer viewer;
    
    int rank;
    MPI_Comm comm;
    
    
   
    A  = self->A;    
    pc = self->pc;
    is = A->is;
   
    maxits  = self->maxits;
    tol     = self->tol*self->tol;
    nloc    = parms_MapGetLocalSize(is);
    comm    = is->comm;

    MPI_Comm_rank(comm, &rank);

    /* --- first permute x and y for pARMS matrix structure ------*/
    parms_VecPerm(x, is); 
    parms_VecPerm(y, is);
/* ---- set isvecperm to true -----*/
    is->isvecperm = true;
    
/*----- allocate memory for working local arrays -------*/
    size = nloc;
    PARMS_NEWARRAY(r0, size);  
    PARMS_NEWARRAY(r,  size);  
    s = r;
    PARMS_NEWARRAY(st, size);
    PARMS_NEWARRAY(t, size);
    PARMS_NEWARRAY0(v, size);
    PARMS_NEWARRAY0(p,  size);
    PARMS_NEWARRAY(pt, size);
    
   
    /*"Knoll trick:"
      parms_PCApply(pc, y, r);*/
   
    /*compute r = A * x */ 
    parms_MatVec(A, x, r);

    /* compute residual vector r = y - Ax */
    parms_VecAYPX(r, y, -1.0, is);    
    PARMS_MEMCPY(p,  r, nloc);
    parms_VecDOT(r, r, &omega, is);

    /* choose r0 (arbitrary) */
    /*PARMS_MEMCPY(r0,r,nloc);*/
    PARMS_MEMCPY(r0, y, nloc);
    parms_VecDOT(r0,r,&rho, is);  
   
    alpha = 0.;
    sigma = rho_alt = 1.0;
 
    its = 0;   
    while(its < maxits  && omega > tol){                 

      parms_PCApply(pc, p, pt);
      parms_MatVec(A, pt, v);
      parms_VecDOT(r0,v,&alpha,is);

      alpha = rho/alpha;
      parms_VecAXPY(x, pt, alpha, is);	
      parms_VecAXPY(s, v, -alpha, is);
      
      parms_VecDOT(s, s, &omega, is);
      if(omega < TINY*TINY){
	its++;
	continue;
      }

      parms_PCApply(pc, s, st);
      parms_MatVec(A, st, t);
      tmp[0] = GDOT(nloc, t, one, s, one);
      tmp[1] = GDOT(nloc, t, one, t, one);
      if(is->isserial == false)
	MPI_Allreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm);
      
      sigma = tmp[0]/tmp[1];
      
      if(ABS_VALUE(sigma) < TINY){
	if(rank == 0)
	  printf("ERROR: sigma = 0\n");
	break;
      }
      
      parms_VecAXPY(x, st, sigma, is);
      parms_VecAXPY(r, t, -sigma, is);

      parms_VecDOT( r, r, &omega, is);	
      parms_VecDOT(r0, r, &rho,   is);
      
      beta = (rho * alpha) / (rho_alt * sigma);
      rho_alt = rho;
      parms_VecAXPY(p, v, -sigma, is);
      parms_VecAYPX(p, r, beta,   is);       

      its++;
    }
    
    if(its == maxits && omega > tol)
      printf("ERROR: no convergence\n");
    
    free(r0);
    free(r);
    free(p);
    free(st);
    free(t);
    free(v);
    
    /* reset isvecperm and do inverse permutation*/
    is->isvecperm = false; // this comes before inverse permutation
    /* permutes x and y */
    parms_VecInvPerm(x, is); 
    parms_VecInvPerm(y, is);

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
