/* This file contains pARMS interface routines to petsc. Note that this is 
 * only compatible with the latest version of petsc (version 3.1). For older 
 * versions (version 3.0) use parms.c_petsc_3.0.
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "petscksp.h"
#include "parms.h"
#include "protos.h"

static   PC_PARMS  *PREC;    /* user-defined preconditioner context */

#undef __FUNCT__  
#define __FUNCT__ "PCCreate_PARMS"
PetscErrorCode PCCreate_PARMS(PC_PARMS **prec)
{
  PC_PARMS *newpc;
  PetscErrorCode       ierr;

  ierr                     = PetscNew(PC_PARMS,&newpc);CHKERRQ(ierr);

  newpc->map		   = 0;
  newpc->A		   = 0;
  newpc->pc		   = 0;
  
  *prec = newpc;

  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCSetUp_PARMS"
PetscErrorCode PCSetUp_PARMS(PC ctx)
{
  Mat pmat;
  PC *precon = &ctx;
  parms_Map   map;
  parms_Mat   A;
  parms_PC    pc;
  char *iluname, *pname;
  int         *maptmp, *mapptr, npro, rank, low, high;
  int         m, n, lsize, *ia, *ja, *ja1, ncols, pos, i;
  int         length, *im;
  double      fill_fact;
  PetscScalar   *aa, *aa1;
  const PetscInt *cols;
  const PetscScalar *values;

/* pARMS options file name */
  char fname[] = "parms_opts";

  MPI_Comm_size(PETSC_COMM_WORLD, &npro);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

/* ------------ Check to make sure compilations for pARMS and PETSC are consistent -----*/
#if defined(DBL_CMPLX)
  if(sizeof(PetscScalar) == sizeof(double))
  {
      fprintf(stderr," %d: Compilation mismatch: Using complex version of pARMS with real version of PETSC !!\n", rank);
      fprintf(stderr," %d: Compile real version of pARMS to use real version of PETSC \n", rank);      
      MPI_Finalize();
      exit(0);
  }
#else
  if(sizeof(PetscScalar) == sizeof(double _Complex))
  {
      fprintf(stderr," %d: Compilation mismatch: Using real version of pARMS with complex version of PETSC !!\n", rank);
      fprintf(stderr," %d: Compile complex version of pARMS to use complex version of PETSC \n", rank);            
      MPI_Finalize();
      exit(0);
  }
#endif

/* Get preconditioning matrix from petsc and setup 
/  pARMS structs
*/
  PCGetOperators(*precon,PETSC_NULL,&pmat,PETSC_NULL);

  MatGetSize(pmat, &n, &m);
  MatGetOwnershipRange(pmat, &low, &high);
  PetscMalloc((npro+1)*sizeof(int), &mapptr);
  PetscMalloc(n*sizeof(int), &maptmp);
  lsize = high - low;

  MPI_Allgather(&lsize, 1, MPI_INT, mapptr+1, 1, MPI_INT, PETSC_COMM_WORLD);
  mapptr[0] = 1;
  for (i = 0; i < npro; i++) {
    mapptr[i+1] += mapptr[i];
  }
  for (i = 0; i < n; i++) {
    maptmp[i] = i+1;
  }

/* create parms map object */
  if(PREC == NULL)
  	parms_MapCreateFromPtr(&map, n, maptmp, mapptr, MPI_COMM_WORLD, 1, NONINTERLACED);
  else /*---- map re-use data structure */
        map = PREC->map;

/* create parms mat object */
  if(PREC == NULL)
  	parms_MatCreate(&A, map);
  else
  {
/*------ Data structure is being re-used ----*/ 
        A = PREC->A;     
        parms_MatReset(A, DIFFERENT_NONZERO_STRUCTURE); 
  }
/* setup and copy csr data structure for parms */
  PetscMalloc((lsize+1)*sizeof(int), &ia);
  ia[0] = 1;
  length = lsize * 8;
  PetscMalloc(length*sizeof(int), &ja);
  PetscMalloc(length*sizeof(PetscScalar), &aa);

  for (i = low; i < high; i++) 
  {
    pos = ia[i-low]-1;
    MatGetRow(pmat, i, &ncols, &cols, &values);
    ia[i-low+1] = ia[i-low] + ncols;

    if (ia[i-low+1] >= length) 
    {
      length += ncols;
      PetscMalloc(length*sizeof(int), &ja1);
      memcpy(ja1, ja, (ia[i-low]-1)*sizeof(int));
      PetscFree(ja);
      ja = ja1;
      PetscMalloc(length*sizeof(PetscScalar), &aa1);
      memcpy(aa1, aa, (ia[i-low]-1)*sizeof(PetscScalar));
      PetscFree(aa);
      aa = aa1;
    }
    memcpy(&ja[pos], cols,  ncols*sizeof(int));
    memcpy(&aa[pos], values,ncols*sizeof(PetscScalar));
    MatRestoreRow(pmat, i, &ncols, &cols, &values);
  }

/* csr info is for local matrix so initialize im[] locally */
  PetscMalloc(lsize*sizeof(int), &im);
  memcpy(im, &maptmp[mapptr[rank]-1], lsize*sizeof(int));

/* convert csr aa vector from petscscalar to double */
  for(i=0; i<ia[lsize]-1; i++){
	aa[i] = (double)aa[i];
        ja[i] = ja[i]+1;
  }

/* Now copy csr matrix to parms_mat object */
  parms_MatSetValues(A, lsize, im, ia, ja, aa, INSERT);

/* free memory */
  PetscFree(maptmp);  
  PetscFree(mapptr); 
  PetscFree(aa); 
  PetscFree(ja);  
  PetscFree(ia);
  PetscFree(im); 

/* setup parms matrix */
  parms_MatSetup(A);

/* Now create parms preconditioner object based on A. 
   Default precon uses BJ with ILU0. Can be changed 
   per user's preference.
*/
  if(PREC == NULL)
  	parms_PCCreate(&pc, A);
  else /*--- re-use pc data structure ---*/
        pc = PREC->pc;
  
/* Set Preconditioner options based on parameters in 
   parms_opts file 
*/
  parms_SetOptions(pc, fname);
  parms_PCSetup(pc);

/* Print some preconditioner stats */
  parms_PCGetName(pc, &pname);
  parms_PCILUGetName(pc, &iluname);
  parms_PCGetRatio(pc, &fill_fact);
  if (rank == 0) {
    printf("The global preconditioner  %2s %s\n", " ", pname);
    printf("The local preconditioner %4s %s\n", " ", iluname);
    printf("The memory usage %12s %-4.2f\n", "=", fill_fact);
    printf("The number of processors %4s %-4d\n", "=", npro);
  }


/* Setup preconditioner object for petsc */
  if(PREC == NULL)
  {
  	PCCreate_PARMS(&PREC);
  	PREC->map = map;
  	PREC->A = A;
  	PREC->pc = pc;  
  }
  return(0);
}

/* prec - precondition object
   b - right-hand-side vector
   x - solution vector on return (preconditioned vector)
*/
#undef __FUNCT__  
#define __FUNCT__ "PCApply_PARMS"
PetscErrorCode PCApply_PARMS(PC dummy,Vec b,Vec x)
{
  PetscErrorCode  ierr;
  PetscScalar     *x1, *b1;
  FLOAT		  *rhs, *rhs1, *sol, *sol1;
  int		  nloc, i, rank;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
/* get petsc vector */
  ierr = VecGetArray(b, &b1);CHKERRQ(ierr);
  ierr = VecGetArray(x, &x1);CHKERRQ(ierr);

/* get size of local system and allocate memory 
   for rhs and solution vectors
*/
  nloc = parms_MapGetLocalSize(PREC->map);
  PetscMalloc((nloc)*sizeof(FLOAT), &rhs);  
  PetscMalloc((nloc)*sizeof(FLOAT), &sol);  
  PetscMalloc((nloc)*sizeof(FLOAT), &rhs1);  
  PetscMalloc((nloc)*sizeof(FLOAT), &sol1);

/* convert from petsscalar to float */
  for(i=0; i<nloc; i++)
  {
    rhs1[i] = (FLOAT)b1[i];
    sol1[i] = (FLOAT)x1[i];  
  }

  parms_VecPermAux(rhs1, rhs, PREC->map);

/* Do preconditioner solve */
  parms_PCApply(PREC->pc, rhs, sol1);

  parms_VecInvPermAux(sol1, sol, PREC->map);

/* convert from float to petscscalar*/
  for(i=0; i<nloc; i++)
  {
    x1[i] = (PetscScalar)sol[i];  
  }

/* restore petsc array and return solution */
  ierr = VecRestoreArray(b, &b1);CHKERRQ(ierr);
  ierr = VecRestoreArray(x, &x1);CHKERRQ(ierr);

  PetscFree(rhs);
  PetscFree(sol);

  return(0);
}

#undef __FUNCT__  
#define __FUNCT__ "PCDestroy_PARMS"
PetscErrorCode PCDestroy_PARMS(PC pc)
{
  PetscErrorCode      ierr;

  parms_MatFree(&(PREC->A));
  parms_MapFree(&(PREC->map));
  parms_PCFree(&(PREC->pc));

  ierr = PetscFree(PREC);CHKERRQ(ierr);
  return(0);
}


/* Some auxilliary functions 
============================
*/

/* Set pARMS preconditioner options. It is easier 
   to simply read from file, so options are set 
   in the file parms_opts
*/
int parms_SetOptions(parms_PC pc, char fname[])
{
  FILE *fp;
  char buf[MAX_LINE], *p, *p1, *p2;
  int i, its, par, lfil[7], nlev, meth[8];
  double eps, droptol[7], tolind;
   
  for(i=0; i<8; i++)
     meth[i] = 0;
 
  if (NULL == (fp = fopen(fname, "r"))) {
    fprintf(stderr, "cannot open parms options file %s\n", fname);
    exit(1);
  }
  
/*----Get parallel Precon Type ---*/
  if(fgets(buf, MAX_LINE, fp)== NULL){
    fprintf(stderr, "Error reading global precon type");
    exit(1);
  }
  STR2U(p, p1);
  if (!strcmp(p1, "RAS")) {
  	parms_PCSetType(pc, PCRAS);
  }
  else if (!strcmp(p1, "SCHUR")){
  	parms_PCSetType(pc, PCSCHUR);
  }
  else if (!strcmp(p1, "BJ")) {
  	parms_PCSetType(pc, PCBJ);
  }
  free(p1);
/* --- Get Local Precon type ---*/
  if(fgets(buf, MAX_LINE, fp) == NULL){
    fprintf(stderr, "Error reading local precon type");
    exit(1);
  }
  STR2U(p, p1);
  if (!strcmp(p1, "ILU0")) {
  	parms_PCSetILUType(pc, PCILU0);
  }
  else if (!strcmp(p1, "ILUK")) {
  	parms_PCSetILUType(pc, PCILUK);
  }
  else if (!strcmp(p1, "ILUT")) {
  	parms_PCSetILUType(pc, PCILUT);
  }
  else if (!strcmp(p1, "ARMS")) {
  	parms_PCSetILUType(pc, PCARMS);
  }
  free(p1);
/*--- Get tolerance for local solve---*/
  if(fscanf(fp,"%lf",&eps) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetInnerEps(pc, eps);
  while(fgetc(fp) != '\n');
/*--- Get number of levels ---*/
  if(fscanf(fp,"%d",&nlev) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetNlevels(pc, nlev);
  while(fgetc(fp) != '\n');  
/*--- Use symmetric or nonsymmetric perm ---*/
  if(fscanf(fp,"%d",&par) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetPermType(pc, par);
  while(fgetc(fp) != '\n');  
/*--- B block size ---*/
  if(fscanf(fp,"%d",&par) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetBsize(pc, par);
  while(fgetc(fp) != '\n');  
/*--- Tolerance for independent sets ---*/
  if(fscanf(fp,"%lf",&tolind) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetTolInd(pc, tolind);
  while(fgetc(fp) != '\n');  
/*--- Inner Krylov dimension ---*/
  if(fscanf(fp,"%d",&par) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetInnerKSize(pc, par);
  while(fgetc(fp) != '\n');  
/*--- maximum number of inner iterations ---*/
  if(fscanf(fp,"%d",&its) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  parms_PCSetInnerMaxits(pc, its);
  while(fgetc(fp) != '\n');
/*--- nonsymmetric perm. for interlevel blocks ---*/
  if(fscanf(fp,"%d",&meth[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');    
/*--- Column Permutations - interlevel blocks ---*/
  if(fscanf(fp,"%d",&meth[1]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n'); 
/*--- Row scaling - interlevel blocks ---*/
  if(fscanf(fp,"%d",&meth[2]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');   
/*--- Column scaling - interlevel blocks ---*/
  if(fscanf(fp,"%d",&meth[3]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');  
/*--- Nonsymm. Perm. for last level SC. ---*/
  if(fscanf(fp,"%d",&meth[4]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');  
/*--- Col. Perm. for last level SC. ---*/
  if(fscanf(fp,"%d",&meth[5]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');  
/*--- Row scaling - last level SC. ---*/
  if(fscanf(fp,"%d",&meth[6]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');  
/*--- Column scaling - last level SC. ---*/
  if(fscanf(fp,"%d",&meth[7]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');  
/*--- lfil0 for ilut, iluk, and arms for lfil[0-3] ---*/
  if(fscanf(fp,"%d",&lfil[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');  
/*--- lfil for shur ---*/
  if(fscanf(fp,"%d",&lfil[4]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');
/*--- lfil for ILUT L and U ---*/
  if(fscanf(fp,"%d",&lfil[5]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');
/*--- droptol0(droptol[0=3], L, U, L^{-1}F, EU^{-1}---*/
  if(fscanf(fp,"%lf",&droptol[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');
/*--- droptol4(for schur complements at each level) ---*/
  if(fscanf(fp,"%lf",&droptol[4]) !=1 ){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');
/*--- droptol5(for ILUT in last level Schur Complement) ---*/
  if(fscanf(fp,"%lf",&droptol[5]) !=1 ){
    fprintf(stderr, "Error reading value");
    exit(1);
  }  
  while(fgetc(fp) != '\n');

  droptol[1] = droptol[2] = droptol[3] = droptol[0]; 
  droptol[6] = droptol[5];
  
  lfil[1] = lfil[2] = lfil[3] = lfil[0];
  lfil[6] = lfil[5];  

  parms_PCSetPermScalOptions(pc, &meth[0], 1);
  parms_PCSetPermScalOptions(pc, &meth[4], 0);

  parms_PCSetFill(pc, lfil);
  parms_PCSetTol(pc, droptol);
  
  fclose(fp);

  return 0;
}
