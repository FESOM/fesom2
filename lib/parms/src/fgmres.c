/*--------------------------------------------------------------------
  fgmres_create    : create the FGMRES solver. 
  fgmres_free      : free the memory for the FGMRES solver. 
  fgmres_setksize  : set the restart size of the FGMRES solver.
  fgmres_view      : dump the FGMRES solver.
  parms_fgmres     : the FGMRES solve function.

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

typedef struct fgmres_data {
  int restart;
  int neigs;
} *fgmres_data;

#define TINY 1.0e-20
int parms_fgmres(parms_Solver self, FLOAT *y, FLOAT *x)
{
    BOOL outflag, intflag;
    int i, i1, pti, pti1, ptih, j, restart, its;
    int ii, jj, k, k1, maxits, nloc, size, one = 1;
    FLOAT  *vv, *z, *hh, *s, *rs, alpha, t1;
    REAL  eps1, tol, ro, t, *c;
#if defined(DBL_CMPLX)
    FLOAT rot;
#else
    REAL gam;    
#endif    

    parms_Mat A;
    parms_PC  pc;
    parms_Map is;
    fgmres_data fdata;
    
    int rank;
    MPI_Comm comm;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    A  = self->A;
    pc = self->pc;
    is = A->is;
    fdata = (fgmres_data)self->data;
    restart = fdata->restart;
    maxits  = self->maxits;
    tol     = self->tol;
    eps1 = tol;
    nloc    = parms_MapGetLocalSize(is);
    comm  =  is->comm;

/* --- first permute x and y for pARMS matrix structure ------*/
    parms_VecPerm(x, is); 
    parms_VecPerm(y, is);
/* ---- set isvecperm to true -----*/
    is->isvecperm = true;

/*----- allocate memory for working local arrays -------*/

/* ---------allocate size for nx(m+1) matrix vv */
    size = nloc*(restart+1); 
    PARMS_NEWARRAY(vv, size);
/*-------- allocate size for nxm matrix z -----*/    
    size = nloc*restart; 
    PARMS_NEWARRAY(z,  size);

/*----- allocate memory for Hessenberg matrix (nonzeros only) 
 *----- and rotation vectors s and rs -------------------*/    
    size = (((restart+1)*(restart+2)/2) - 1) + 2*(restart+1);
    PARMS_NEWARRAY(hh, size);

    s = hh + ((restart+1)*(restart+2)/2) - 1;
    rs = s + restart + 1;
/*--------- allocate memory for rotation matrix c -------
 * This is done separately since for complex valued systems
 * c is still a real-valued vector, and hence cannot be 
 * allocated as sizeof(complex double) as is the case for 
 * s and rs and hh above --------------------------------*/    
    PARMS_NEWARRAY(c, restart+1);
    
    /* outer loop starts here */
    its = 0;

    outflag = true;
    while(outflag) {
      /* compute vv[0] = A * x */
      parms_MatVec(A, x, vv);

      /* compute residual vector */
      parms_VecAYPX(vv, y, -1.0, is);
      /* compute the norm of the residual */
      parms_VecGetNorm2(vv, &ro, is);

      if(ABS_VALUE(ro) <= DBL_EPSILON)  {      
	outflag = false;
	break;
      }
      t1 = 1.0 / ro;
      parms_VecScale(vv, t1, is);      
      if(its == 0)
	eps1 = tol*ro;

/* ----------initialize 1-st term of rhs of hessenberg system ------------*/

      rs[0] = ro;

      i = -1;
      pti = 0;
      pti1 = 0;
      ptih = 0;
      intflag = true;
      while (intflag) {
	i++;
	its++;
	i1 = i + 1;
        pti = i*nloc;
        pti1 = i1*nloc;
        
/*------------- preconditioning operation z = K^{-1}vv ---------------*/
	parms_PCApply(pc, &vv[pti], &z[pti]);

/*------------- compute A*z -----------------*/
	parms_MatVec(A, &z[pti], &vv[pti1]);
	
/*------------- classical Gram - Schmidt -------------------*/
#if defined(DBL_CMPLX)
/* ----- Check for serial case ------ */
        if(is->isserial)
        {
            for(j=0; j<i1; j++)
            {
	      hh[ptih+j]  = 0.;
	      for (i=0;i<nloc;i++) hh[ptih+j] += vv[i+j*nloc] * vv[pti1+i];

                 alpha = -hh[ptih+j];
		 for (i=0; i<nloc; i++) vv[pti1+i] += alpha * vv[j*nloc+i];
      	    }
      	 }
      	 else /*-------------- parallel case -------*/
      	 {
            for(j=0; j<i1; j++)
	      {
		t1 = 0.;
                for (i=0;i<nloc;i++) t1 += vv[i+j*nloc] * vv[pti1+i];
                MPI_Allreduce(&t1, &hh[ptih+j], one, MPI_CMPLX, MPI_CMPLX_SUM, comm);
                alpha = -hh[ptih+j];
		for (i=0; i<nloc; i++) vv[pti1+i] += alpha * vv[j*nloc+i];

      	    }
      	}       	           	 
#else
/* ----- Check for serial case ------ */
        if(is->isserial)
        {
            for(j=0; j<i1; j++)
            {
              hh[ptih+j]  = 0.;
              for (i=0;i<nloc;i++) hh[ptih+j] += vv[i+j*nloc] * vv[pti1+i];
                 alpha = -hh[ptih+j];
                 for (i=0; i<nloc; i++) vv[pti1+i] += alpha * vv[j*nloc+i];
      	    }
      	 }
      	 else /*-------------- parallel case -------*/
      	 {
            for(j=0; j<i1; j++)
            {
              t1 = 0.;
	      for (i=0;i<nloc;i++) t1 += vv[i+j*nloc] * vv[pti1+i];

                MPI_Allreduce(&t1, &hh[ptih+j], one, MPI_DOUBLE, MPI_SUM, comm);
                alpha = -hh[ptih+j];
		for (i=0; i<nloc; i++) vv[pti1+i] += alpha * vv[j*nloc+i];

      	    }
      	}  
#endif  
        parms_VecGetNorm2(&vv[pti1], &t, is);      
        hh[ptih+i1] = t;

	if (fabs(t) > TINY) {
	  t1 = 1.0 / t;
	  parms_VecScale(&vv[pti1], t1, is);
	}

	/* done with classical Gram-Schmidt and Arnoldi step. now update
	 * factorization of hh */
#if defined(DBL_CMPLX)		 
	if (i != 0) {
	  for(k = 1; k <= i; k++) {
	    k1 = k-1;
	    t1 = hh[ptih+k1];
	    
	    hh[ptih+k1] = c[k1]*t1 + s[k1]*hh[ptih+k];
	    hh[ptih+k] = -conj(s[k1])*t1 + c[k1]*hh[ptih+k];	    
	  }
      }
/*-----------get next plane rotation------------ */
      zclartg(hh[ptih+i], hh[ptih+i1], &c[i], &s[i], &rot);
      rs[i1] = -conj(s[i])*rs[i];
      rs[i] = c[i]*rs[i];      
      hh[ptih+i] = rot;
      ro = cabs(rs[i1]);	
#else
	if (i != 0) {
	  for(k=1; k<=i; k++){ 
	    k1 = k-1;
	    t1 = hh[ptih+k1];
	    
	    hh[ptih+k1] = c[k1]*t1 + s[k1]*hh[ptih+k];
	    hh[ptih+k] = -s[k1]*t1 + c[k1]*hh[ptih+k];
	  }
	}
/*-----------get next plane rotation------------ */
	gam = sqrt(hh[ptih+i]*hh[ptih+i] + hh[ptih+i1]*hh[ptih+i1]);
      /* 
	 if gamma is zero then any small value will do ...
	 will affect only residual estimate 
      */	
	if (fabs(gam) <= TINY) gam = TINY;
	
	/* determine-next-plane-rotation */
	c[i] = hh[ptih+i]/gam;
	s[i] = hh[ptih+i1]/gam;
	
	rs[i1] = -s[i]*rs[i];
	rs[i] = c[i]*rs[i];
	/* determine res. norm and test for convergence */
	hh[ptih+i] = c[i]*hh[ptih+i] + s[i]*hh[ptih+i1];
	ro = fabs(rs[i1]);	
#endif	

/*------------ Check for convergence ---------*/    
	if ((i+1 >= restart) || (ro <= eps1) || its >= maxits)   
	   intflag = false;
	else
/*------------ update hh pointer ptih ---------*/   
           ptih += i+2;    	
      }

      /* now compute solution first solve upper triangular system */
      rs[i] = rs[i]/hh[ptih+i];
      for (ii = 1; ii <= i; ii++) {
	k = i-ii;
	k1 = k+1;
	t1 = rs[k];
	for (j = k1; j <= i; j++) {
	  jj = ((j+1)*(j+2)/2) - 1;
	  t1 = t1 - hh[jj+k]*rs[j];
	}
	jj = ((k+1)*(k+2)/2)-1;
	rs[k] = t1/hh[jj+k];
      }
	/* done with back substitution. now form linear combination to
	 * get solution */
      for (j = 0; j <= i; j++) {
	t1 = rs[j];
	parms_VecAXPY(x, &z[j*nloc], t1, is);
      }
      /* test for return */
      /*      if ((ro <= eps1) || (its >= maxits)) (AF){ */
      if((ro <= tol) || (its >= maxits)) {
	outflag = false;
      }
    }

    free(vv);
    free(z);
    free(hh);
    free(c);


/* reset isvecperm and do inverse permutation*/
    is->isvecperm = false; // this comes before inverse permutation
    /* permutes x and y */
    parms_VecInvPerm(x, is); 
    parms_VecInvPerm(y, is);

    self->its = its;
    return 0;
}

static int fgmres_getresidual(parms_Solver self, FLOAT *y, FLOAT *x, FLOAT *res)
{

    parms_Mat A;

    A  = self->A;

    parms_MatMVPY(A, -1.0, x, 1.0, y, res);

    return 0;
}
     
static int fgmres_getresidualnorm2(parms_Solver self, FLOAT *y, FLOAT *x, REAL *rnorm)
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


static int fgmres_free(parms_Solver *self)
{
  fgmres_data fdata;

/* free data associated with solver */  
  fdata = (fgmres_data)(*self)->data;
  PARMS_FREE(fdata);  
  return 0;
}

static int fgmres_view(parms_Solver self, parms_Viewer v)
{
  FILE *fp;
  char *name, *iluname;
  int restart;
  fgmres_data fdata;

  fdata = (fgmres_data)self->data;
  restart = fdata->restart;

  parms_ViewerGetFP(v, &fp);

  fprintf(fp,"\n=============================\n");
  fprintf(fp,"	Solver Parameters	\n");
  fprintf(fp,"=============================\n");
  
  fprintf(fp, "Solver type = flexible gmres (fgmres) \n");
  
  fprintf(fp, "maxits = %d \n", self->maxits);
  fprintf(fp, "Relative tolerance = %-8.2e \n", self->tol);

  fprintf(fp, "Krylov dimension = %d \n", restart);

  parms_PCGetName(self->pc, &name);
  parms_PCILUGetName(self->pc, &iluname);
  fprintf(fp, "Global Preconditioner: %s\n", name);
  fprintf(fp, "Local Preconditioner: %s\n", iluname);

  parms_ViewerGetFP(v, &fp);

  return 0;
}

static int fgmres_setksize(parms_Solver self, int restart)
{
  fgmres_data fdata;

  fdata = (fgmres_data)self->data;
  fdata->restart = restart;

  return 0;
}

static int fgmres_setneig(parms_Solver self, int neigs)
{
  fgmres_data fdata;

  fdata = (fgmres_data)self->data;
  fdata->neigs = neigs;

  return 0;
}

static struct parms_Solver_ops parms_fgmres_sol = {
  parms_fgmres,
  fgmres_getresidual,
  fgmres_getresidualnorm2,
  fgmres_setksize,
  fgmres_setneig,
  fgmres_free,
  fgmres_view
};


/** Create the FGMRES solver. 
 *
 *  \param self A parms_Solver object.
 *  \return 0 on success.
 */
int fgmres_create(parms_Solver self)
{
  fgmres_data fdata;

  PARMS_NEW(fdata);
  fdata->restart = 30;
  fdata->neigs = 0;
  self->data =fdata;
  PARMS_MEMCPY(self->ops, &parms_fgmres_sol, 1);		

  return 0;
}

