/*--------------------------------------------------------------------
  cg_create    : create the CG solver. 
  cg_free      : free the memory for the CG solver. 
  cg_view      : dump the CG solver.
  parms_cg     : the CG solve function.

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

typedef struct cg_data {
  int neigs;
} *cg_data;

#define TINY 1.0e-20
int parms_cg(parms_Solver self, FLOAT *y, FLOAT *x)
{
  int i, its, j;
    int maxits, nloc;
    FLOAT alpha, tmp[2];
    FLOAT *r, *z, *p, *Ap;
    REAL  tol, omega, beta, r2;
   
    
    parms_Mat A;
    parms_PC  pc;
    parms_Map is;
    parms_vcsr Acsr;
    
    MPI_Comm comm;
    
       
    A  = self->A;    
    pc = self->pc;
    is = A->is;
   
    maxits  = self->maxits;
    tol     = self->tol * self->tol;
    nloc    = parms_MapGetLocalSize(is);
    comm    = is->comm;
    

/*----- allocate memory for local working arrays -------*/

    r =  malloc(nloc*sizeof(FLOAT));  
    z =  malloc(nloc*sizeof(FLOAT));  
    p =  malloc(nloc*sizeof(FLOAT));
    Ap = malloc(nloc*sizeof(FLOAT));
    
/* --- first permute x and y for pARMS matrix structure ------*/
    if (is->isperm) { 
        for (i = 0; i < nloc; i++) {
	  z[is->perm[i]] = x[i];
	  r[is->perm[i]] = y[i];
	}
	memcpy(x, z, nloc*sizeof(FLOAT));
	memcpy(y, r, nloc*sizeof(FLOAT));
	is->isvecperm = true;
    }

    omega = 2.*tol;

    /* compute residual vector r = y - Ax */
    parms_MatVec(A,x,r);
    for (i = 0; i < nloc; i++)  r[i] = y[i] - r[i];
    
    parms_PCApply(pc,r,z);

    r2 = 0;
    for (i = 0; i < nloc; i++){
      p[i] = z[i];
      r2 += r[i] * z[i];
    }

    if(is->isserial == false)
      MPI_Allreduce(MPI_IN_PLACE, &r2, 1, MPI_DOUBLE, MPI_SUM, comm);
    

    for (its = 1; its <= maxits; its++){                 
      
       parms_MatVec(A,p,Ap);

      alpha = 0.;      
      for (i = 0; i < nloc; i++) alpha += p[i] * Ap[i];            
      if(is->isserial == false)
	MPI_Allreduce(MPI_IN_PLACE, &alpha, 1, MPI_DOUBLE, MPI_SUM, comm);

      alpha = r2 / alpha;
      
      for (i = 0; i < nloc; i++){
	x[i] += alpha *  p[i];
	r[i] -= alpha * Ap[i];
      }

      parms_PCApply(pc, r, z);

      tmp[0] = 0.;
      tmp[1] = 0.;
      for (i = 0; i < nloc; i++){
	tmp[0] += r[i] * r[i];  /* omega */
	tmp[1] += r[i] * z[i];
      }
      if(is->isserial == false)
	MPI_Allreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm);

      omega = tmp[0];
      /* Convergence criterion reached? */
      if (omega <= tol) break;  

      beta = tmp[1] / r2;
      r2   = tmp[1];
      
      for (i = 0; i < nloc; i++){
	p[i] = z[i] + beta *  p[i];
      }
    }
    
    if(its == maxits && omega > tol)
      printf("ERROR: no convergence\n");
    
    /* reset isvecperm and do inverse permutation*/
    if (is->isperm) { 
        for (i = 0; i < nloc; i++) {
	  r[is->iperm[i]] = x[i];
	  z[is->iperm[i]] = y[i];
	}
	memcpy(x, r, nloc*sizeof(FLOAT));
        memcpy(y, z, nloc*sizeof(FLOAT));
	is->isvecperm = false; 
    }

    free(z);  
    free(r);  
    free(p);
    free(Ap);


    self->its = its;
    return 0;
}

static int cg_getresidual(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *res)
{

    parms_Mat A;

    A  = self->A;

    parms_MatMVPY(A, -1.0, x, 1.0, y, res);

    return 0;
}
     
static int cg_getresidualnorm2(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm)
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


static int cg_free(parms_Solver *self)
{
  /*dummy*/
  return 0;
}

static int cg_view(parms_Solver self, parms_Viewer v)
{
  FILE *fp;
  char *name, *iluname;
  int dim;
  parms_ViewerGetFP(v, &fp);

  fprintf(fp,"\n=============================\n");
  fprintf(fp,"	Solver Parameters	\n");
  fprintf(fp,"=============================\n");
  
  fprintf(fp, "Solver type = Conjugate Gradient Method (CG) \n");
  
  fprintf(fp, "maxits = %d \n", self->maxits);
  fprintf(fp, "Relative tolerance = %-8.2e \n", self->tol);

  parms_PCGetName(self->pc, &name);
  parms_PCILUGetName(self->pc, &iluname);
  fprintf(fp, "Global Preconditioner: %s\n", name);
  fprintf(fp, "Local Preconditioner: %s\n", iluname);

  parms_ViewerGetFP(v, &fp);

  return 0;
}

static int cg_setksize(parms_Solver self, int restart)
{
  /*dummy*/
  return 0;
}

static int cg_setneig(parms_Solver self, int neigs)
{
  /*dummy*/
  return 0;
}

static struct parms_Solver_ops parms_cg_sol = {
  parms_cg,
  cg_getresidual,
  cg_getresidualnorm2,
  cg_setksize,
  cg_setneig,
  cg_free,
  cg_view
};


/** Create the CG solver. 
 *
 *  \param self A parms_Solver object.
 *  \return 0 on success.
 */
int cg_create(parms_Solver self)
{
  PARMS_MEMCPY(self->ops, &parms_cg_sol, 1);		
  return 0;
}
