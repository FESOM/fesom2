/*-----------------------------------------------------------------------
 * Test program
 WORKS ONLY WHEN  nprocx * nprocy == np 
 *-----------------------------------------------------------------------
 *
 *  Each node generates  its own  part of the  5-point  matrix using a
 *  simple regular partitioning.  Then the rhs  is generated, the local
 *  data structure is set and the system is  solved using
 *  preconditioned FGMRES.
 *
 *-----------------------------------------------------------------------
 *
 * the mesh is  defined as follows. The grid is virtually laid out on a
 * 2-dimensional array of processors: mprocx in x direction, mprocy  in
 * y direction.  mprocx is  read from  a  file  (inputs)  and mprocy is
 * determined as mprocy = numproc / mprocx  where  numproc is the total
 * number of processors available for this run. The simplest case is to
 * always take mprocx=1.  Then the mesh sizeso in each direction are set
 * by the commands:
 * 
 *      nx = nmesh*mprocx
 *      ny = nmesh*mprocy
 *      nz = 1
 *---------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "parms.h"
#include "aux.h"

int main(int argc, char *argv[])
{
  /* n      -- number of total variables
   * nloc   -- number of local variables 
   * nx     -- number of variables in the x direction
   * ny     -- number of variables in the y direction
   * nz     -- 1
   * iov    -- overlap (always 0 for vertex based partitioning)
   * mprocx -- number of processors in the x direction
   * mprocy -- number of processors in the y direction
   * maxits -- maximum number of iterations
   * nmesh  -- mesh points in x and y directions
   * im     -- the dimension of Krylov subspace
   * lfil   -- lfil for ilut
   * maptmp -- global labelings which are stored processor by processor
   * mapptr -- a pointer array. mapptr[i] points to the beginning of
   *           processor i
   * nnz    -- the number of none zero elements
   */
  int npro, myid,n,nloc,nx,ny,nz,iov,mprocx,mprocy,mprocz,xnmesh,ynmesh, nnz; 
  int *maptmp,*mapptr, iout, *ja, *ia, *lvars, size,nrow,start;
  int its;
  double tpc, ttol, res0, res1, ratio;
  int count;
  fprm prm;
  char *name, *iluname;
  FILE *fp;
  /* variables related to pARMS */
  parms_Map    map;
  FLOAT    *sol, *rhs, *y;
  parms_Mat    A;
  parms_PC     pc;
  parms_Solver solver;
  parms_Timer  t;

/* Viewer object  */
//  parms_Viewer  mv; 

  /* 
   * a      -- a in CSR format
   * stencil-- working array
   */
  FLOAT stencil[100], *a;

/*-------------------- Initialize MPI environment */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npro);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

/*--------- read data from file input ------------*/
  prm = malloc(sizeof(*prm));
  if (prm == NULL) {
    fprintf(stderr, "cannot allocate memory for prm\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }

  if (argc >= 2) {
    read_param(argv[1], prm);
  }
  else if (argc == 1){
    read_param("inputs", prm);
  }

  if (myid == 0) {
    fprintf(stderr, "mprocx = %d, mprocy = %d\n", prm->mprocx, prm->mprocy);
    fprintf(stderr, "xnmesh = %d, ynmesh = %d\n", prm->xnmesh, prm->ynmesh);
  }

  mprocx = prm->mprocx;
  mprocy = prm->mprocy;
  xnmesh = prm->xnmesh;
  ynmesh = prm->ynmesh;

/*------------------- mprocy = dm->comm->npro / mprocx; */
  if (npro !=  mprocx*mprocy) {
    if (myid == 0) {
  fprintf(stderr, " ** ERROR nproc is not equal to mprocx*mprocy = %d\n",
	  mprocx*mprocy);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 99);
    exit(1);
  }

  nx = mprocx * xnmesh;
  ny = mprocy * ynmesh;
  nz = 1;
  size = xnmesh*ynmesh*mprocx*mprocy;

  mapptr = malloc((npro+1)*sizeof(*mapptr));
  if (mapptr == NULL) {
    fprintf(stderr, "cannot allocate memory for mapptr\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }
  maptmp = malloc(size*sizeof(*maptmp));
  if (maptmp == NULL) {
    fprintf(stderr, "cannot allocate memory for maptmp\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }

  iov = 0;
  mprocz = 1;
/*-------------------- call part2_ from skit.f to generate local mesh */
  part2_(&nx, &ny, &nz, &mprocx, &mprocy, &mprocz, &iov, maptmp, 
	mapptr, &iout);
  n = nx * ny * nz;
/*-------------------- Create map object and setup mapping from 
                       global labels to  processors */
  parms_MapCreateFromPtr(&map, n, maptmp, mapptr, MPI_COMM_WORLD, 1, NONINTERLACED);
  nloc = parms_MapGetLocalSize(map);
  lvars = (int *)malloc(nloc*sizeof(int));
  parms_MapGetGlobalIndices(map, lvars);
/*-------------------- Allocate temporary arrays for storing 
                       local matrix in CSR format   */ 
  nnz = 10*xnmesh*ynmesh;

  a = malloc(nnz*sizeof(*a));
  if (a == NULL) {
    fprintf(stderr, "cannot allocate memory for a\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }
  ja = malloc(nnz*sizeof(*ja));
  if (ja == NULL) {
    fprintf(stderr, "cannot allocate memory for ja\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }
  ia = malloc((nloc+1)*sizeof(*ia));
  if (ia == NULL) {
    fprintf(stderr, "cannot allocate memory for ia\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }

  nrow = mapptr[myid+1] - mapptr[myid];
  start = mapptr[myid] - mapptr[0];
/*-------------------- Generate the local matrix using gen5loc (in fdmat.f) */
#if defined(DBL_CMPLX)
  zgen5loc_(&nx, &ny, &nz, &nloc, &maptmp[start],a,ja,ia,stencil); 
#else
  gen5loc_(&nx, &ny, &nz, &nloc, &maptmp[start],a,ja,ia,stencil); 
#endif  
/*-------------------- Create a matrix object based on map */
  parms_MatCreate(&A, map);
  
/*-------------------- Insert entries into the matrix */
  parms_MatSetValues(A, nrow, &maptmp[start], ia, ja, a, INSERT); 
/*-------------------- Free the matrix stored in CSR format (ia,ja,a) */
  free(a);
  free(ja);
  free(ia);
  
  free(mapptr);
  free(maptmp);  
/*-------------------- setup matrix A */
  parms_MatSetup(A);
/*-------------------- create timer */
  parms_TimerCreate(&t);
/*-------------------- create preconditioner */
  parms_PCCreate(&pc, A);
/*-------------------- set parameters for pc */
  set_pc_params(pc, prm);
/*-------------------- reset the timer */
  parms_TimerReset(t);
/*-------------------- setup preconditioner */
  parms_PCSetup(pc);
  tpc = parms_TimerGet(t);
/*-------------------- stats: get ratio of nnz(preconditioning matrix) 
                       over nnz(A) */
  parms_PCGetRatio(pc, &ratio);
/*-------------------- pause the timer */
  parms_TimerPause(t);
/*-------------------- Create a krylov subspace object */
  parms_SolverCreate(&solver, A, pc);
/*-------------------- Set pamaters for solver */
  set_solver_params(solver, prm);

  free(prm);
/*-------------------- Create solution and right-hand-side vectors */

  rhs = (FLOAT *)malloc(nloc*sizeof(FLOAT));
  sol = (FLOAT *)malloc(nloc*sizeof(FLOAT));
  y = (FLOAT *)malloc(nloc*sizeof(FLOAT));
  
/*------------- print matrix --- DEBUG --------*/
/*
  parms_ViewerCreate(&mv, "matrix.out");
  parms_MatViewCOO_dvcsr(A, mv);  
*/

/*-------------------- Set up an artifical right-hand-side vector */  
  for(count=0; count<nloc; count++)
  	sol[count] = 1.0;

  parms_MatVec(A, sol, rhs);
/*-------------------- Get the initial guess */

  for(count=0; count<nloc; count++)
  	sol[count] = 0.0;
  
  parms_MatVec(A, sol, y);
  parms_VecAXPY(y, rhs, -1.0, map);
/*-------------------- Get the initial residual norm */
  parms_VecGetNorm2(y, &res0, map);
/*-------------------- solve linear systems with the krylov subspace 
                       method selected */
  parms_TimerRestart(t);

  parms_SolverApply(solver, rhs, sol);
/*-------------------- total elapsed time */
  ttol = parms_TimerGet(t);
/*--------------------  number iterations */
  its = parms_SolverGetIts(solver);
/*-------------------- Compute the residual error */
  parms_SolverGetResidualNorm2(solver, rhs, sol, &res1);

/*-------------------- processor 0 outputs the result */
  if (myid == 0) {
    if(argc == 3) {
      if((fp = fopen(argv[2], "w")) == NULL) {
	fprintf(stderr, "cannot open file %s\n", argv[2]);
	exit(1);
      }
    }
    else {
      fp = stdout;
    }
    char buf[40];
    parms_PCGetName(pc, &name);
    parms_PCILUGetName(pc, &iluname);
    fprintf(fp, "The preconditioner  %9s %s\n", " ", name);
    fprintf(fp, "The local preconditioner %4s %s\n", " ", iluname);
    fprintf(fp, "The memory usage %12s %-4.2f\n", "=", ratio);
    fprintf(fp, "The number of processors %4s %-4d\n", "=", npro);
    fprintf(fp, "The number of iterations %4s %-4d\n", "=", its);
    sprintf(buf, "%8.2fs", tpc);
    fprintf(fp, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
    sprintf(buf, "%8.2fs", ttol-tpc);
    fprintf(fp, "The solving time %12s %-s\n", "=", strtok(buf, " "));
    sprintf(buf, "%8.2fs", ttol);
    fprintf(fp, "The total time %14s %-s\n", "=", strtok(buf, " "));
    fprintf(fp, "The initial residual error  = %-8.2e\n", res0);
    fprintf(fp, "The solution residual error = %-8.2e\n", res1);

    if(fp != stdout) {
      fclose(fp);
    }
  }


/*-------------------- free memories for map, vec, mat, solver */


  free(sol);
  free(rhs);
  free(y);
  parms_MatFree(&A);
  parms_MapFree(&map);
  free(lvars);
  parms_PCFree(&pc);
  parms_SolverFree(&solver);
  parms_TimerFree(&t);

/*  
  parms_ViewerFree(&mv);
*/
/*--------------------  exit MPI environment */
  MPI_Finalize();
  return 0;
}
