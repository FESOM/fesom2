#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "aux.h"

int read_param(char *fname, char mname[MAX_MAT][MAX_LINE], fprm prm)
{
  FILE *fp, *fmat;
  char buf[MAX_LINE], *tmpstring;
  char *p, *p1, *p2;
  int i; 

  if (NULL == (fp = fopen(fname, "r"))) {
    fprintf(stderr, "cannot open file %s\n", fname);
    exit(1);
  }

  for (i = 0; i < 18; i++) {
    prm->ipar[i] = -1;
  }
  
/*----Get parallel Precon Type ---*/
  if(fgets(buf, MAX_LINE, fp) == NULL){
    fprintf(stderr, "Error reading global precon type");
    exit(1);
  }
  STR2U(p, p1);
  if (!strcmp(p1, "RAS")) {
  	prm->pctype = PCRAS;
  }
  else if (!strcmp(p1, "SCHUR")){
	prm->pctype = PCSCHUR;
  }
  else if (!strcmp(p1, "BJ")) {
  	prm->pctype = PCBJ;
  }
  free(p1);
/* --- Get Local Precon type ---*/
  if(fgets(buf, MAX_LINE, fp) == NULL){
    fprintf(stderr, "Error reading local precon type");
    exit(1);
  }
  STR2U(p, p1);
  if (!strcmp(p1, "ILU0")) {
	prm->pcilutype = PCILU0;
  }
  else if (!strcmp(p1, "ILUK")) {
	prm->pcilutype = PCILUK;
  }
  else if (!strcmp(p1, "ILUT")) {
	prm->pcilutype = PCILUT;
  }
  else if (!strcmp(p1, "ARMS")) {
	prm->pcilutype = PCARMS;
  }
  free(p1);
/*--- Get tolerance for local solve---*/
  if(fscanf(fp,"%lf",&prm->pgfpar[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- Get tolerance for global solve ---*/
  if(fscanf(fp,"%lf",&prm->pgfpar[1]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- Get number of levels ---*/
  if(fscanf(fp,"%d",&prm->ipar[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Use symmetric or nonsymmetric perm ---*/
  if(fscanf(fp,"%d",&prm->ipar[1]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- B block size ---*/
  if(fscanf(fp,"%d",&prm->ipar[2]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Tolerance for independent sets ---*/
  if(fscanf(fp,"%lf",&prm->tolind) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- Outer Krylov dimension ---*/
  if(fscanf(fp,"%d",&prm->ipar[6]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');    
/*--- Maximum outer iterations ---*/
  if(fscanf(fp,"%d",&prm->ipar[7]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Inner Krylov dimension ---*/
  if(fscanf(fp,"%d",&prm->ipar[4]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- maximum number of inner iterations ---*/
  if(fscanf(fp,"%d",&prm->ipar[5]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- nonsymmetric perm. for interlevel blocks ---*/
  if(fscanf(fp,"%d",&prm->ipar[10]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');    
/*--- Column Permutations - interlevel blocks ---*/
  if(fscanf(fp,"%d",&prm->ipar[11]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n'); 
/*--- Row scaling - interlevel blocks ---*/
  if(fscanf(fp,"%d",&prm->ipar[12]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');   
/*--- Column scaling - interlevel blocks ---*/
  if(fscanf(fp,"%d",&prm->ipar[13]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Nonsymm. Perm. for last level SC. ---*/
  if(fscanf(fp,"%d",&prm->ipar[14]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Col. Perm. for last level SC. ---*/
  if(fscanf(fp,"%d",&prm->ipar[15]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Row scaling - last level SC. ---*/
  if(fscanf(fp,"%d",&prm->ipar[16]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- Column scaling - last level SC. ---*/
  if(fscanf(fp,"%d",&prm->ipar[17]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- lfil0 for ilut, iluk, and arms for lfil[0-3] ---*/
  if(fscanf(fp,"%d",&prm->lfil[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');  
/*--- lfil for shur ---*/
  if(fscanf(fp,"%d",&prm->lfil[4]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- lfil for ILUT L and U ---*/
  if(fscanf(fp,"%d",&prm->lfil[5]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- droptol0(droptol[0=3], L, U, L^{-1}F, EU^{-1}---*/
  if(fscanf(fp,"%lf",&prm->droptol[0]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- droptol4(for schur complements at each level) ---*/
  if(fscanf(fp,"%lf",&prm->droptol[4]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    
  while(fgetc(fp) != '\n');
/*--- droptol5(for ILUT in last level Schur Complement) ---*/
  if(fscanf(fp,"%lf",&prm->droptol[5]) != 1){
    fprintf(stderr, "Error reading value");
    exit(1);
  }    

  prm->droptol[1] = prm->droptol[2] = prm->droptol[3] =
    prm->droptol[0]; 
  prm->droptol[6] = prm->droptol[5];
  
  prm->lfil[1] = prm->lfil[2] = prm->lfil[3] =
    prm->lfil[0];
  prm->lfil[6] = prm->lfil[5];  
  
  fclose(fp);
    
/*------- Done reading inputs, now read list of Matrices -------------*/
#if defined(DBL_CMPLX)
  if (NULL == (fmat = fopen("./matfileCmplx", "r"))) {
    fprintf(stderr, "cannot open file: matfile\n");
    exit(1);
  }
#else  
    if (NULL == (fmat = fopen("./matfileReal", "r"))) {
    fprintf(stderr, "cannot open file: matfile\n");
    exit(1);
  }
#endif  
   i = 0;
    do {
       if(fgets(buf,MAX_LINE,fmat) == NULL){
         fprintf(stderr, "Error reading matrix list");
         exit(1);
       }
       tmpstring = (char *) strtok(buf, " ");
       tmpstring++;  
       tmpstring[strlen(tmpstring)-1]='\0';
      strcpy(mname[i], --tmpstring);
      i++;
    } while (strncmp(buf, "##", 2));

    fclose(fmat);
 
  return 0;
}
  
void set_pc_params(parms_PC pc, fprm prm)
{
  parms_PCSetType(pc,            prm->pctype);
  parms_PCSetILUType(pc,         prm->pcilutype);
  parms_PCSetNlevels(pc,         prm->ipar[0]);
  parms_PCSetPermType(pc,        prm->ipar[1]);
  parms_PCSetBsize(pc,           prm->ipar[2]);
  parms_PCSetInnerEps(pc,        prm->pgfpar[0]);
  parms_PCSetInnerKSize(pc,      prm->ipar[4]);
  parms_PCSetInnerMaxits(pc,     prm->ipar[5]);
  parms_PCSetFill(pc,            prm->lfil);
  parms_PCSetTol(pc,             prm->droptol);
  parms_PCSetTolInd(pc,          prm->tolind);
  parms_PCSetPermScalOptions(pc, &prm->ipar[10], 1);
  parms_PCSetPermScalOptions(pc, &prm->ipar[14], 0);
}
void set_solver_params(parms_Solver solver, fprm prm)
{
  char buf[BUFFLEN];

  sprintf(buf, "%d", prm->ipar[7]);
  parms_SolverSetParam(solver, MAXITS, buf);
  sprintf(buf, "%d", prm->ipar[6]);
  parms_SolverSetParam(solver, KSIZE,  buf);
  sprintf(buf, "%g", prm->pgfpar[1]);
  parms_SolverSetParam(solver, DTOL,  buf);
}

void fread_param_(char *fname, fprm *prm, char *matrix, int *matlen, int len)
{
  char mname[MAX_MAT][MAX_LINE], *buff, *buff2;

  buff2 = malloc((len+1)*sizeof(*buff2));
  strncpy(buff2, fname, len);
  buff2[len] = '\0';
  *prm = malloc(sizeof(**prm));
  read_param(buff2, mname, *prm);
  buff = &mname[0][0];
  strncpy(matrix, buff, strlen(buff));
  matrix[strlen(buff)] = '\0';
  *matlen = (int)strlen(matrix);
  free(buff2);

}

void fset_pc_params_(parms_PC *pc, fprm *prm)
{
  set_pc_params(*pc, *prm);
}

void fset_solver_params_(parms_Solver *solver, fprm *prm)
{
  set_solver_params(*solver, *prm);
}
  
void fprm_free_(fprm *prm)
{
  free(*prm);
}
