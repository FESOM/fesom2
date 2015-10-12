/*      Program driver - preconditioners in pARMS  can be called
 *      through PETSc.
 *----------------------------------------------------------------------
 *  In this test program the matrix is read from file.
 *  The matrix is assumed to be in Harwell-Boeing
 *  format. The problem is solved with PETSc. The preconditioners in
 *  pARMS can be used through PETSc.
 **/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "petscao.h"
#include "petscksp.h"
#include "parms.h"
#include "protos.h"

int main(int argc, char *argv[])
{
  int myid,usmy,ierr,len;
  /* declarations related to Harwell-boeing format for reading the HB
     matri. Second part is related to I/O parameters */
  char matrix[BUFLEN], guesol[2], title[72],
    key[8], type[3]; 
  int nrhs, nc, n, nnz, tmp0, tmp, tmp2, tmp3, job;

  /* working array for reading matrix */
  FLOAT *a,*rhstmp,*b;
  int *ja, *ia,*jb,*ib;
  int npro, *iwk, i1,i;
  int mnnz, ii, jj;

  /* working array for symmetryzing matrix */
  FLOAT *mc;
  int *jc, *ic;
  int nnzrow, pos;

  Vec       sol, rhs, y;
  Mat       A;
  KSP       ksp;
  PC        pc;
  PetscReal res1, res2, value;
  int       j, iters;

/* User defined preconditioner */
  PetscTruth     user_defined_pc;  /* flag option for user defined preconditioner */

  /* initialize MPI environment */
  PetscInitialize(&argc, &argv, (char *)0, (char *)0);
  MPI_Comm_rank(PETSC_COMM_WORLD, &myid);
  MPI_Comm_size(PETSC_COMM_WORLD, &npro);

/* ------------ Check to make sure compilations for pARMS and PETSC are consistent -----*/
#if defined(DBL_CMPLX)
  if(sizeof(PetscScalar) == sizeof(double))
  {
      fprintf(stderr," %d: Compilation mismatch: Using complex version of pARMS with real version of PETSC !!\n", myid);
      fprintf(stderr," %d: Compile real version of pARMS to use real version of PETSC \n", myid);      
      MPI_Finalize();
      exit(0);
  }
#else
  if(sizeof(PetscScalar) == sizeof(double _Complex))
  {
      fprintf(stderr," %d: Compilation mismatch: Using real version of pARMS with complex version of PETSC !!\n", myid);
      fprintf(stderr," %d: Compile complex version of pARMS to use complex version of PETSC \n", myid);            
      MPI_Finalize();
      exit(0);
  }
#endif

  tmp0 = 0;
  nrhs = 1;
  /* Read matrix; either using user-defined function (routine uread) or 
     |               SPARSKIT function for reading Harwell-Boeieng matrices
  */
  /* read matrix name from file input */
  read_matrix(matrix);
  
  /* variable "matrix" stores the name of the file in HB format 
     |
     |   Read a Harwell-Boeing matrix. using wreadmtc c-version of
     |      sparsit routine - call wreadmtc a first time to determine sizes
     |      of arryas. read in values on the second call. 
  */
  len = strlen(matrix);
  a = NULL; ja = NULL; ia = NULL; rhstmp = NULL;
#if defined(DBL_CMPLX)  
  zwreadmtc_(&tmp0,&tmp0,&tmp0,matrix,&len,a,ja,ia,rhstmp,&nrhs,
	   guesol,&n,&nc,&nnz,title,key,type,&ierr); 
#else
  wreadmtc_(&tmp0,&tmp0,&tmp0,matrix,&len,a,ja,ia,rhstmp,&nrhs,
	   guesol,&n,&nc,&nnz,title,key,type,&ierr); 
#endif  
  ia = malloc((n+1)*sizeof(*ia));
  ja = malloc(nnz*sizeof(*ja));
  a  = malloc(nnz*sizeof(*a));
  
  if (nrhs) {
    tmp = 3;
    tmp2 = n;
    tmp3 = nnz;
    if (guesol[0] == 'G' || guesol[0] == 'g') nrhs++;
    if (guesol[1] == 'X' || guesol[1] == 'x') nrhs++;
    rhstmp = malloc((n*nrhs)*sizeof(*rhstmp));
  }
  if (!nrhs) {
    tmp = 2;
    tmp2 = n;
    tmp3 = nnz;
  }
  
  /* Array sizes determined. Now call wreadmtc again for really
     reading */
#if defined(DBL_CMPLX)     
  zwreadmtc_(&tmp2,&tmp3,&tmp,matrix,&len,a,ja,ia,rhstmp,&nrhs,
	   guesol,&n,&nc,&nnz,title,key,type,&ierr);
#else
  wreadmtc_(&tmp2,&tmp3,&tmp,matrix,&len,a,ja,ia,rhstmp,&nrhs,
	   guesol,&n,&nc,&nnz,title,key,type,&ierr);
#endif	   
  free(rhstmp);
  if(ierr > 0 && ierr < 4 && !myid) {
    fprintf(stderr, "cannot read matrix\n");
    fprintf(stderr, "ierr = %d\n", ierr);
    exit(1);
  }
  else if (ierr >= 4 && !myid){
    fprintf(stderr, "*** Warning: RHS or solution vectors are not read. Artificial RHS will be used\n");
    nrhs = 0;
  }
  
  if(myid == 0) {
    fprintf(stdout,"READ the matrix %.*s %.*s \n",8,key,3,type);
  }
  

  if (myid == 0){
    fprintf(stdout,"Matrix dimension is %d, Number of nonzeros is %d\n",n,nnz);
  }
  
  
  /* the matrix is symmetrized */
  b  = malloc(nnz*sizeof(*b));
  jb = malloc(nnz*sizeof(*jb));
  ib = malloc((n+1)*sizeof(*ib));
  
  job = 1;
  i = 1;
/* ---- convert from CSC to CSR matrix ------*/  
#if defined(DBL_CMPLX)  
  zcsrcsc_(&n, &job, &i, a, ja, ia, b, jb, ib);
#else
  csrcsc_(&n, &job, &i, a, ja, ia, b, jb, ib);
#endif

/*----- symmetrize or not ---- */
  usmy = 0;
  if(usmy == 1) {
    iwk = malloc(n*sizeof(*iwk));
    mc  = malloc((2*nnz)*sizeof(*mc));
    jc  = malloc((2*nnz)*sizeof(*jc));
    ic  = malloc((n+1)*sizeof(*ic));
    
    for(i = 0; i < nnz; i++) {
      b[i] = 0.0;
    }
    job = 1; 
    i1 = 2*nnz;
    
    /* compute C = A + B */
#if defined(DBL_CMPLX)    
    zaplb_(&n, &n, &job, a, ja, ia, b, jb, ib, mc, jc, ic, &i1, iwk,
	 &ierr);
#else
    aplb_(&n, &n, &job, a, ja, ia, b, jb, ib, mc, jc, ic, &i1, iwk,
	 &ierr);
#endif	 
    if (ierr) fprintf(stderr, "Error: in aplb, ierr = %d\n", ierr);
    nnz = ic[n]-1;
    if (!myid)
      fprintf(stderr,
	      "Matrix pattern has been symmetrized; added %d nonzeros.\n",
	      ic[n]-ia[n]);

    ja = realloc(ja, nnz*sizeof(*ja));
    a  = realloc(a,  nnz*sizeof(*a));

    memcpy(ia, ic, (n+1)*sizeof(*ia));
    memcpy(ja, jc, nnz*sizeof(*ia));
    memcpy(a, mc, nnz*sizeof(*a));

    free(b);
    free(jb);
    free(ib);
    free(mc);
    free(jc);
    free(ic);
    free(iwk);
  }
  else
  {
    memcpy(ia, ib, (n+1)*sizeof(*ia));
    memcpy(ja, jb, nnz*sizeof(*ia));
    memcpy(a, b, nnz*sizeof(*a));  
    free(b);
    free(jb);
    free(ib);  
  }

  /* create Mat object */
  mnnz = 0;
  for (i = 1; i <= n; i++) {
    nnzrow = ia[i] - ia[i-1];
    if (mnnz < nnzrow) {
      mnnz = nnzrow;
    }
  }

  /* create a matrix object */
  MatCreateMPIAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, n, n, mnnz, PETSC_NULL, mnnz,
		  PETSC_NULL, &A);

  for(i = 1; i <= n; i++) {
    nnzrow = ia[i] - ia[i-1];
    pos = ia[i-1] - ia[0];
    ii = i-1;
    for (j = 0; j < nnzrow; j++) {
      jj = ja[pos+j] - 1;
      MatSetValues(A, 1, &ii, 1, &jj, &a[pos+j], INSERT_VALUES);
    }
  }

  /* free temporary arrays */
  free(a); free(ja); free(ia);

  /* assemble the matrix */
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  /* create an abstract krylov subspace object */
  KSPCreate(PETSC_COMM_WORLD, &ksp);

  /* create solution and right-hand-side vector */
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &rhs);
  VecDuplicate(rhs, &sol);
  VecDuplicate(rhs, &y);

  value = 1.0;
  VecSet(sol, value);
  MatMult(A, sol, rhs); 

  /* set initial guess */
  value = 0.0;
  VecSet(sol, value);

  VecAssemblyBegin(sol);
  VecAssemblyEnd(sol);

  /* initial residual norm */
  MatMult(A, sol, y);
  value = -1.0;
  VecAXPY(rhs, value, y);
  VecNorm(y, NORM_2, &res1);

  /* solve linear systems with the krylov subspace method selected */
  KSPSetOperators(ksp, A, A, SAME_NONZERO_PATTERN);
  KSPSetType(ksp, KSPFGMRES); 
  KSPGMRESSetRestart(ksp, 20);
  KSPSetTolerances(ksp, 1.0e-6, 1.0e-20, 1.0e6, 1000);
  /* Extract PC context from the KSP context */
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  /*
     Set a user-defined "shell" preconditioner if desired
  */
  ierr = PetscOptionsHasName(PETSC_NULL,"-user_defined_pc",&user_defined_pc);CHKERRQ(ierr);
  if (user_defined_pc) {

    /* (Required) Indicate to PETSc that we're using a "shell" preconditioner */
    ierr = PCSetType(pc,PCSHELL);CHKERRQ(ierr);

    /* (Required) Set the user-defined routine for applying the preconditioner */
    ierr = PCShellSetApply(pc,PCApply_PARMS);CHKERRQ(ierr);
    ierr = PCShellSetContext(pc,&pc);CHKERRQ(ierr);

    /* (Optional) Set a name for the preconditioner, used for PCView() */
    ierr = PCShellSetName(pc,"pARMS Preconditioner");CHKERRQ(ierr);

    /* (Optional) Do any setup required for the preconditioner */
    ierr = PCShellSetSetUp(pc,PCSetUp_PARMS);CHKERRQ(ierr);

    /* (Optional) Destructor for the preconditioner */
    ierr = PCShellSetDestroy(pc,PCDestroy_PARMS);CHKERRQ(ierr);

  } else {
    ierr = PCSetType(pc,PCBJACOBI);CHKERRQ(ierr);
  }

  KSPSetFromOptions(ksp);

  KSPSolve(ksp, rhs, sol);
  KSPGetIterationNumber(ksp, &iters);
  VecNorm(rhs, NORM_2, &res1);
  KSPGetResidualNorm(ksp, &res2);
  if(myid == 0) {
    fprintf(stderr, "The number of iteration %5s %d\n", "=",iters);
    fprintf(stderr, "The initial residual error  = %8.3e\n", res1);
    fprintf(stderr, "The solution residual error = %8.3e\n", res2);
  }

  /* free memories for map, vec, mat, solver */
  MatDestroy(A);
  VecDestroy(sol);
  VecDestroy(rhs);
  VecDestroy(y);
  KSPDestroy(ksp);

  /* exit MPI environment */
  PetscFinalize();
  return EXIT_SUCCESS;
}
