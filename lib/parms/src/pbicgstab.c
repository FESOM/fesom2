/*--------------------------------------------------------------------
  pbicgstab_create    : create the PBICGS solver. 
  pbicgstab_free      : free the memory for the PBICGS solver. 
  pbicgstab_view      : dump the PBICGS solver.
  parms_pbicgstab     : the PBICGS solve function.

  Pipelined BiCGstab, overlapping global sums with computation. See:
  Siegfried Cools, Wim Vanroose
  The communication-hiding pipelined BiCGStab method for the parallel solution of large unsymmetric linear systems
  Parallel Computing 65, pp. 1-20, July 2017
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

typedef struct pbicgs_data {
  int neigs;
} *pbicgstab_data;

#define TINY 1.0e-20
int parms_pbicgstab(parms_Solver self, double *rhs, double *xin)
{

    int i, its;
    int maxits, nloc, size, one = 1;
    double tmp[5];
    double *b, *r0, *r, *rt, *w, *wt, *t, *pt, *s, *st, *z, *zt;
    double *q, *qt, *v, *x, *y;
    double  tol, alpha, beta, omega, myres, rho, rho_alt;
    MPI_Request req, req1, req2;

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
/*----- initialize to zero if not initialized before the iteration starts ---*/
    size = nloc;

    PARMS_NEWARRAY(b, size); 
    PARMS_NEWARRAY(r0, size); 
    PARMS_NEWARRAY(r, size);  for (i=0; i<nloc; i++) r[i] = 0.;
    PARMS_NEWARRAY(rt, size); 
    PARMS_NEWARRAY(w, size);  
    PARMS_NEWARRAY(wt, size); 
    PARMS_NEWARRAY(t, size);  
    PARMS_NEWARRAY(pt, size); for (i=0; i<nloc; i++) pt[i] = 0.;
    PARMS_NEWARRAY(s, size);  for (i=0; i<nloc; i++) s[i] = 0.;
    PARMS_NEWARRAY(st, size); for (i=0; i<nloc; i++) st[i] = 0.;
    PARMS_NEWARRAY(z, size);  for (i=0; i<nloc; i++) z[i] = 0.;
    PARMS_NEWARRAY(zt, size); for (i=0; i<nloc; i++) zt[i] = 0.;
    PARMS_NEWARRAY(q, size);  for (i=0; i<nloc; i++) q[i] = 0.;
    PARMS_NEWARRAY(qt, size); for (i=0; i<nloc; i++) qt[i] = 0.;
    PARMS_NEWARRAY(y, size);  for (i=0; i<nloc; i++) y[i] = 0.;
    PARMS_NEWARRAY(v, size);  for (i=0; i<nloc; i++) v[i] = 0.;
    PARMS_NEWARRAY(x, size);

/* --- first permute x and y for pARMS matrix structure ------*/
    if (is->isperm) { 
      for (i = 0; i < nloc; i++) {
	x[is->perm[i]] = xin[i];
	b[is->perm[i]] = rhs[i];
      }
      is->isvecperm = true;
    }
   
    
    /* compute residual vector r0 = b - Ax */
    /* b = y  (permuted), t=Ax temporary vector */
    parms_MatVec(A, x, t);
    for (i=0; i<nloc; i++) r0[i] = b[i] - t[i];  
    
    parms_PCApply(pc, r0, rt);  // rt = M r0
    parms_MatVec(A, rt, w);      // w = A rt

    tmp[0] = 0.; tmp[1] = 0.;
    for (i=0; i<nloc; i++) {
      tmp[0] += r0[i]*r0[i];
      tmp[1] += r0[i]*w[i];
        r[i]  = r0[i];
    }
    if(is->isserial == false)
      MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req); 

    parms_PCApply(pc, w, wt);   // wt = M w
    parms_MatVec(A, wt, t);     //  t = Awt

    omega = 0.;
    beta = 0.;

    if(is->isserial == false) MPI_Wait(&req, MPI_STATUS_IGNORE);

    rho_alt = tmp[0];
    alpha = tmp[0] / tmp[1];  // alpha = (r0,r0) / (r0,w)


    for (its = 0; its < maxits; its++){                 

      for (i=0; i<nloc; i++) pt[i] = rt[i] + beta*(pt[i] - omega*st[i]);
      for (i=0; i<nloc; i++)  s[i] =  w[i] + beta*( s[i] - omega* z[i]);
      for (i=0; i<nloc; i++) st[i] = wt[i] + beta*(st[i] - omega*zt[i]);
      for (i=0; i<nloc; i++)  z[i] =  t[i] + beta*( z[i] - omega* v[i]);

      for (i=0; i<nloc; i++)  q[i] =  r[i] - alpha* s[i];
      for (i=0; i<nloc; i++) qt[i] = rt[i] - alpha*st[i];
      for (i=0; i<nloc; i++)  y[i] =  w[i] - alpha* z[i];

      tmp[0] = 0.; tmp[1] = 0.;
      for (i=0; i<nloc; i++) {
	tmp[0] += q[i]*y[i];
	tmp[1] += y[i]*y[i];
      }
      if(is->isserial == false)
	MPI_Iallreduce(MPI_IN_PLACE, tmp, 2, MPI_DOUBLE, MPI_SUM, comm, &req1); 
      
      parms_PCApply(pc, z, zt);
      parms_MatVec(A, zt, v);  

      if(is->isserial == false) MPI_Wait(&req1, MPI_STATUS_IGNORE);
      
      omega = tmp[0]/tmp[1];

      for (i=0; i<nloc; i++)  x[i] += alpha*pt[i] + omega*qt[i];
      for (i=0; i<nloc; i++)  r[i]  =  q[i] - omega*y[i];
      for (i=0; i<nloc; i++) rt[i]  = qt[i] - omega*(wt[i] - alpha*zt[i]);
      for (i=0; i<nloc; i++)  w[i]  =  y[i] - omega*( t[i] - alpha* v[i]);

      tmp[0] = 0.; tmp[1] = 0.; tmp[2] = 0.; tmp[3] = 0.; tmp[4] = 0.;
      for (i=0; i<nloc; i++) {
	tmp[0] += r0[i]*r[i];
	tmp[1] += r0[i]*w[i];
	tmp[2] += r0[i]*s[i];
	tmp[3] += r0[i]*z[i];
	tmp[4] +=  r[i]*r[i]; // residuum
      }
      if(is->isserial == false)
	MPI_Iallreduce(MPI_IN_PLACE, tmp, 5, MPI_DOUBLE, MPI_SUM, comm, &req2);
            
      parms_PCApply(pc, w, wt);
      parms_MatVec(A, wt, t);  
 
      if(is->isserial == false) MPI_Wait(&req2, MPI_STATUS_IGNORE);
      
      // Convergence check:
      if ( tmp[4] < tol) break;

      rho = tmp[0];
      beta = alpha*rho / (omega*rho_alt);
      rho_alt = rho;
      alpha = rho / (tmp[1] + beta*(tmp[2] - omega*tmp[3])); 

      //pbicgstab_getresidualnorm2(self, b, x, &myres);
      
      //if(rank == 0) printf("nparms check at %d, myres=%e, (r,r)=%e\n",its,myres, tmp[4]);
      
      //if (myres < tol) break;
    }
    

    its++;
    if(rank == 0){
      if(its == maxits) {
	printf("ERROR: no convergence\n");
      } else {
	printf("Parms iterations, residuum %d, %e\n",its,sqrt(tmp[4]));
      }
    }
    /* reset isvecperm and do inverse permutation*/
    if (is->isperm) { 
        for (i = 0; i < nloc; i++) {
	  xin[is->iperm[i]] = x[i];
	}
	is->isvecperm = false; 
    }

    
    free(b);
    free(r0);
    free(r);
    free(rt);
    free(w);
    free(wt);
    free(t);
    free(pt);
    free(s);
    free(st);
    free(z);
    free(zt);
    free(q);
    free(qt);
    free(y);
    free(v);
    free(x);

    self->its = its;
    return 0;
}

static int pbicgstab_getresidual(parms_Solver self, double *y, double *x, double *res)
{

    parms_Mat A;

    A  = self->A;

    parms_MatMVPY(A, -1.0, x, 1.0, y, res);

    return 0;
}
     
static int pbicgstab_getresidualnorm2(parms_Solver self, double *y, double *x, double *rnorm)
{

    int nloc;
    double *res;
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


static int pbicgstab_free(parms_Solver *self)
{
  /*dummy*/
  return 0;
}

static int pbicgstab_view(parms_Solver self, parms_Viewer v)
{
  FILE *fp;
  char *name, *iluname;
  int dim;
  parms_ViewerGetFP(v, &fp);

  fprintf(fp,"\n=============================\n");
  fprintf(fp,"	Solver Parameters	\n");
  fprintf(fp,"=============================\n");
  
  fprintf(fp, "Solver type = Pipelined BiConjugate Gradient Stabilized Method (PBICGS) \n");
  
  fprintf(fp, "maxits = %d \n", self->maxits);
  fprintf(fp, "Relative tolerance = %-8.2e \n", self->tol);

  parms_PCGetName(self->pc, &name);
  parms_PCILUGetName(self->pc, &iluname);
  fprintf(fp, "Global Preconditioner: %s\n", name);
  fprintf(fp, "Local Preconditioner: %s\n", iluname);

  parms_ViewerGetFP(v, &fp);

  return 0;
}

static int pbicgstab_setksize(parms_Solver self, int restart)
{
  /*dummy*/
  return 0;
}

static int pbicgstab_setneig(parms_Solver self, int neigs)
{
  /*dummy*/
  return 0;
}

static struct parms_Solver_ops parms_pbicgstab_sol = {
  parms_pbicgstab,
  pbicgstab_getresidual,
  pbicgstab_getresidualnorm2,
  pbicgstab_setksize,
  pbicgstab_setneig,
  pbicgstab_free,
  pbicgstab_view
};


/** Create the PBICGS solver. 
 *
 *  \param self A parms_Solver object.
 *  \return 0 on success.
 */
int pbicgstab_create(parms_Solver self)
{
  PARMS_MEMCPY(self->ops, &parms_pbicgstab_sol, 1);		
  return 0;
}
