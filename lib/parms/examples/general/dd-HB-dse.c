/*---------------------------------------------------------------------- 
 *                           Program dd-HB-dse  
 *----------------------------------------------------------------------
 *
 *  In this test program, each processor reads the whole matrix
 *  from file. The matrix is assumed to be in Harwell-Boeing format.
 *  Matrix graph is then  partitioned  using  DSE, a simple partitioning
 *  routine, and scatters the local matrices to each processor. Once
 *  these submatrices are received each processor solves the problem
 *  using preconditioned FGMRES preconditioned with :
 *                       BJ, RAS, SCHUR
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(__ICC)
#include <mathimf.h>
#else
#include <math.h>
#endif 
#include "parms.h"
#include "aux.h"

int main(int argc, char *argv[])
{

  /* declarations related to Harwell-boeing format for reading the HB 
     matri. Second part is related to I/O parameters */
  char mname[MAX_MAT][MAX_LINE], guesol[2], title[72], key[8], type[3];
  int nrhs, nc, n, nnz, tmp0, tmp, tmp2, tmp3, job, mat;
  int myid, ierr, i, nloc;
  /* memory usage of the preconditioning matrix */
  double ratio;
  /* working array for reading matrix */
  double norm, res1, tpc, ttol;
  FLOAT *a, *rhstmp;
  int *ja, *ia;
  int npro,its, *im;
  fprm prm;
  FILE *fp=NULL, *fout=NULL;
  char *name, *iluname, buf[40];

/*-------------------- variables related to dse partitioning */
  int *riord, *dom, *idom, *mask, *jwk, *link;

/*-------------------- variables related to pARMS */
  parms_Map       map;
  FLOAT       *x, *y, *rhs, *resvec;
  parms_Mat       A;
  parms_PC        pc;
  parms_Solver    solver;
  parms_Timer     tm;
 
/* Viewer object for solver */
//  parms_Viewer  sv;  

/*-------------------- initialize MPI environment */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &npro);

  tmp0 = 0;
  nrhs = 0;

/*-------------------- read matrix name from input file */
  prm = malloc(sizeof(*prm));
  if (prm == NULL) {
    fprintf(stderr, "cannot allocate memory for prm\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }

  if (argc >= 2) {
    read_param(argv[1],mname, prm);
  }
  else if (argc == 1){
    read_param("inputs",mname, prm);
  }

  /* variable "mname" stores the name of the file in HB format. Read a
     Harwell-Boeing matrix. using wreadmtc c-version of sparsit
     routine - call wreadmtc a first time to determine sizes of
     arrys. read in values on the second call.  
  */

/* --- Begin loop over matrices ----*/
  mat = 0;
  while(mname[mat][1] != '#'){
  a = NULL; ja = NULL; ia = NULL; rhstmp = NULL;

#if defined(DBL_CMPLX)
  zreadmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
	    guesol,&n,&nc,&nnz,title,key,type,&ierr);  
#else
  readmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
	    guesol,&n,&nc,&nnz,title,key,type,&ierr);   
#endif

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
  ia = malloc((n+1)*sizeof(*ia));
  if (ia == NULL) {
    fprintf(stderr, "cannot allocate memory for ia\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }
  rhstmp = malloc(n*sizeof(*rhstmp));
  if (rhstmp == NULL) {
    fprintf(stderr, "cannot allocate memory for rhstmp\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }
  if(nrhs != 0)
    tmp = 3;
  else
    tmp = 2;
  tmp2 = n;
  tmp3 = nnz;

  /*-------------------- Array sizes determined. Now call 
                         wreadmtc again for really reading */
/*
  wreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],&len,a,ja,ia,rhstmp,&nrhs,
	    guesol,&n,&nc,&nnz,title,key,type,&ierr);
*/
#if defined(DBL_CMPLX)
  zreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
	    guesol,&n,&nc,&nnz,title,key,type,&ierr);
#else
  readmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
	    guesol,&n,&nc,&nnz,title,key,type,&ierr);
#endif	    

  if(ierr != 0) {
    fprintf(stderr, "ierr = %d\n", ierr);
    fprintf(stderr, "cannot read matrix\n");
    MPI_Finalize();
    exit(1);
  }

  if(myid == 0){
    if(argc == 3) {
      if (NULL == (fp = fopen(argv[2], "w"))) {
	fprintf(stderr, "cannot open file %s\n", argv[2]);
	exit(1);
      }
    }
    else {
      fp = stdout;
    }
    fprintf(fp, "\nMatrix %d: %.*s %.*s \n",(mat+1),8,key,3,type);
    fprintf(fp, "n = %d, nnz = %d\n", n, nnz);
  }
  
/*-------------Convert from CSC to CSR format ------------*/
    int *jb, *ib;
    FLOAT *b;
    b   = malloc(nnz*sizeof(*b));
    if (b == NULL) {
      fprintf(stderr, "cannot allocate memory for b\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    jb  = malloc(nnz*sizeof(*jb));
    if (jb == NULL) {
      fprintf(stderr, "cannot allocate memory for jb\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    ib  = malloc((n+1)*sizeof(*ib));
    if (ib == NULL) {
      fprintf(stderr, "cannot allocate memory for ib\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    job = 1;
#if defined(DBL_CMPLX)
    zcsrcsc_(&n, &job, &job, a, ja, ia, b, jb, ib);
#else
    csrcsc_(&n, &job, &job, a, ja, ia, b, jb, ib);
#endif
/*---------------Copy CSR matrix ------------------*/

    memcpy(ia, ib, (n+1)*sizeof(*ia));
    memcpy(ja, jb, nnz*sizeof(*ja));
    memcpy(a,  b, nnz*sizeof(*a));

/*---------------Free CSR matrix ------------------*/

   free(ib);
   free(b);
   free(jb);  

  idom = malloc((npro+1)*sizeof(*idom));
  if (idom == NULL) {
    fprintf(stderr, "cannot allocate memory for idom\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }
  dom = malloc(n*sizeof(*dom));
  if (dom == NULL) {
    fprintf(stderr, "cannot allocate memory for dom\n");
    MPI_Abort(MPI_COMM_WORLD, 66);
  }

  if (npro == 1) {
    for (i = 0; i < n; i++) {
      dom[i] = i+1;
    }
    idom[0] = 1;
    idom[1] = n+1;
  }
  else {
    riord = malloc(n*sizeof(*riord));
    if (riord == NULL) {
      fprintf(stderr, "cannot allocate memory for riord\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    mask = malloc(n*sizeof(*mask));
    if (mask == NULL) {
      fprintf(stderr, "cannot allocate memory for mask\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    jwk = malloc(2*n*sizeof(*jwk));
    if (jwk == NULL) {
      fprintf(stderr, "cannot allocate memory for jwk\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    link = malloc(n*sizeof(*link));
    if (link == NULL) {
      fprintf(stderr, "cannot allocate memory for link\n");
      MPI_Abort(MPI_COMM_WORLD, 66);
    }
    dse_(&n, ja, ia, &npro, riord, dom, idom, mask, jwk, link);
    free(riord);
    free(mask);
    free(jwk);
    free(link);
  }
/*
    for (i = 0; i < n; i++) {
      dom[i] = i+1;
    }
    for(i=1; i<npro; i++)
	idom[i] = idom[i-1]+(n/npro); 
    idom[npro+1] = n+1;
printf("idom[%d] = %d; dom[%d] = %d \n", myid, idom[myid], idom[myid]-1,dom[idom[myid]-1]);
*/
/*-------------------- Create map object */
  parms_MapCreateFromPtr(&map, n, dom, idom, MPI_COMM_WORLD, 1, NONINTERLACED);

  nloc = parms_MapGetLocalSize(map);
  
  /* Free dom and idom */
  free(dom);
  free(idom);

/*-------------------- Create distributed vectors based on map */  
  x = (FLOAT *)malloc(nloc*sizeof(FLOAT));
  rhs = (FLOAT *)malloc(nloc*sizeof(FLOAT));
  y = (FLOAT *)malloc(nloc*sizeof(FLOAT));
  resvec = (FLOAT *)malloc(nloc*sizeof(FLOAT));

/*-------------------- create a distributed matrix based on map */
  parms_MatCreate(&A, map);
/*-------------------- Insert values into A */
/* Initialize global indices array */
  im = (int *)malloc(n*sizeof(int));
  for (i = 0; i < n; i++) {
    im[i] = i+1;
  }
  parms_MatSetValues(A, n, im, ia, ja, a, INSERT);
//  parms_MatSetElementMatrix(A, n, im, ia, ja, a, INSERT);
//  parms_MatAssembleElementMatrix(A);
/*-------------------- free the matrix stored in CSR format(a,ja,ia) */
  free(a);
  free(ja);
  free(ia);

/*-------------------- Setup the matrix and communication structure */
  parms_MatSetup(A);
/*-------------------- Copy rhs or Setup artifical right-hand-side vector */
//nrhs = 0;
  if(!nrhs)
  {
     for(i=0; i<nloc; i++)
     {
  	x[i] = 1.0;
     }
     parms_MatVec(A, x, rhs);
  }
  else
  {
    parms_VecSetValues(rhs, n, im, rhstmp, INSERT, map);
  }

 free(rhstmp);
 free(im);
/*--------------------Setup initial guess to a vector of all-0. */
  for(i=0; i<nloc; i++)
  {
  	x[i] = 0.0;
  }
/*--------------------Get 2-norm of initial residual */
  parms_MatVec(A, x, resvec);

  for(i=0; i<nloc; i++)
  {
  	resvec[i] = resvec[i] - rhs[i];
  }
  parms_VecGetNorm2(resvec, &norm, map);
  free(resvec);
/*--------------------Create preconditioner based on the matrix A. */
  parms_PCCreate(&pc, A);

/*--------------------set parameters for pc */
  set_pc_params(pc, prm);

/*--------------------create a timer */
  parms_TimerCreate(&tm);
/*--------------------reset the timer */
  parms_TimerReset(tm);
  parms_PCSetup(pc);

/*--------------------get the elapsed time spent on creating PC */
  tpc =  parms_TimerGet(tm);

/*--------------------get the ratio of the number of nonzero entries in the
     preconditioning matrix to that in the original matrix */
  parms_PCGetRatio(pc, &ratio);
/*--------------------pause the timer */
  parms_TimerPause(tm);

/*--------------------Create a solver based on A and pc */

  parms_SolverCreate(&solver, A, pc);
/*--------------------Set the solver type */
  parms_SolverSetType(solver, SOLFGMRES);

/*--------------------set parameters for solver */
  set_solver_params(solver, prm);

/*--------------------set up solver -- no longer needed - DOK */
//  parms_SolverSetup(solver);

/*--------------------restart the timer */
  parms_TimerRestart(tm);

/*--------------------Solve the linear equation */
  parms_SolverApply(solver, rhs, x);

/*--------------------get total time spent on creating the pc and solving the linear
     system */
  ttol = parms_TimerGet(tm);

/*--------------------Get the number of iterations */
  its = parms_SolverGetIts(solver);

/*--------------------Compute the residual error  */
    parms_MatVec(A, x, y);
    for(i=0; i<nloc; i++)
    {
	y[i] = rhs[i] - y[i];
    }
  parms_VecGetNorm2(y, &res1, map);
/*--------------------processor 0 outputs the result */
  if (myid == 0) {
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
    fprintf(fp, "The initial residual norm %3s %-8.2e\n", "=",norm);
    fprintf(fp, "The final residual norm %5s %-8.2e\n","=", res1);
    
    if(fp != stdout)
	fclose(fp);
/* --- Write output to file ---- */
    fout = fopen("output.txt", "aw");
    fprintf(fout, "\nMatrix %d: %.*s %.*s \n",(mat+1),8,key,3,type);
    fprintf(fout, "n = %d, nnz = %d\n", n, nnz);
    fprintf(fout, "The preconditioner  %9s %s\n", " ", name);
    fprintf(fout, "The local preconditioner %4s %s\n", " ", iluname);
    fprintf(fout, "The memory usage %12s %-4.2f\n", "=", ratio);
    fprintf(fout, "The number of processors %4s %-4d\n", "=", npro);
    fprintf(fout, "The number of iterations %4s %-4d\n", "=", its);
    sprintf(buf, "%8.2fs", tpc);
    fprintf(fout, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
    sprintf(buf, "%8.2fs", ttol-tpc);
    fprintf(fout, "The solving time %12s %-s\n", "=", strtok(buf, " "));
    sprintf(buf, "%8.2fs", ttol);
    fprintf(fout, "The total time %14s %-s\n", "=", strtok(buf, " "));
    fprintf(fout, "The initial residual norm %3s %-8.2e\n", "=",norm);
    fprintf(fout, "The final residual norm %5s %-8.2e\n","=", res1);
    fclose(fout);
  }

/*
  parms_ViewerCreate(&sv, "solver.out");
  parms_TimerView(tm, sv);
*/
/*--------------------Free memories */
  free(x);
  free(y);
  free(rhs);
  parms_MatFree(&A);
  parms_MapFree(&map);
  parms_PCFree(&pc);
  parms_SolverFree(&solver);
  parms_TimerFree(&tm);

/*  
  parms_ViewerFree(&sv);
*/
/*----Goto next matrix ---*/
  mat++;

 }

/*--------------Free prm--------------*/
  free(prm);
/*--------------------Exit MPI environment */
  MPI_Finalize();  
  return 0;
}


