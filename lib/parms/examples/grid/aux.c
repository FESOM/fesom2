#include <stdio.h>
#include <ctype.h>
#include "aux.h"

int read_param(char *fname, fprm prm)
{
  FILE *fp;
  char buf[MAXLINE];
  char *p, *p1, *p2;
  int its;

  if (NULL == (fp = fopen(fname, "r"))) {
    fprintf(stderr, "cannot open file %s\n", fname);
    exit(1);
  }
  
  its = 0;
  while (fgets(buf, MAXLINE, fp)) {
    switch (its) {
    case 0:
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
      break;
    case 1:
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
      break;
    case 2:
      sscanf(buf,"%lf",&prm->pgfpar[0]);
      break;
    case 3:
      sscanf(buf,"%lf",&prm->pgfpar[1]);
      break;
    case 4:
      sscanf(buf,"%d",&prm->ipar[0]);
      break;
    case 5:
      sscanf(buf,"%d",&prm->ipar[1]);
      break;
    case 6:
      sscanf(buf, "%d", &prm->ipar[2]);
      break;
    case 7:
      sscanf(buf,"%lf",&prm->tolind);
      break;
    case 8:
      sscanf(buf,"%d",&prm->ipar[6]);
      break;
    case 9:
      sscanf(buf,"%d", &prm->ipar[7]);
      break;
    case 10:
      sscanf(buf, "%d", &prm->ipar[4]);
      break;
    case 11:
      sscanf(buf, "%d", &prm->ipar[5]);
      break;
    case 12:
      sscanf(buf, "%d", &prm->ipar[10]);
      break;
    case 13:
      sscanf(buf, "%d", &prm->ipar[11]);
      break;
    case 14:
      sscanf(buf, "%d", &prm->ipar[12]);
      break;
    case 15:
      sscanf(buf, "%d", &prm->ipar[13]);
      break;
    case 16:
      sscanf(buf, "%d", &prm->ipar[14]);
      break;
    case 17:
      sscanf(buf, "%d", &prm->ipar[15]);
      break;
    case 18:
      sscanf(buf, "%d", &prm->ipar[16]);
      break;
    case 19:
      sscanf(buf, "%d", &prm->ipar[17]);
      break;
    case 20:
      sscanf(buf, "%d", &prm->lfil[0]);
      break;
    case 21:
      sscanf(buf, "%d", &prm->lfil[4]);
      break;
    case 22:
      sscanf(buf, "%d", &prm->lfil[5]);
      break;
    case 23:
      sscanf(buf, "%lf", &prm->droptol[0]);
      break;
    case 24:
      sscanf(buf, "%lf", &prm->droptol[4]);
      break;
    case 25:
      sscanf(buf, "%lf", &prm->droptol[5]);
      break;
    case 26:
      sscanf(buf, "%d", &prm->mprocx);
      break;
    case 27:
      sscanf(buf, "%d", &prm->mprocy);
      break;
    case 28:
      sscanf(buf, "%d", &prm->xnmesh);
      break; 
    case 29:
      sscanf(buf, "%d", &prm->ynmesh);
      break; 
    default:
      break;
    }
    its++;
  }
  fclose(fp);

  prm->droptol[1] = prm->droptol[2] = prm->droptol[3] =
    prm->droptol[0]; 
  prm->droptol[6] = prm->droptol[5];
  
  prm->lfil[1] = prm->lfil[2] = prm->lfil[3] =
    prm->lfil[0];
  prm->lfil[6] = prm->lfil[5];
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
  char buf[MAXLINE];

  sprintf(buf, "%d", prm->ipar[7]);
  parms_SolverSetParam(solver, MAXITS, buf);
  sprintf(buf, "%d", prm->ipar[6]);
  parms_SolverSetParam(solver, KSIZE,  buf);
  sprintf(buf, "%g", prm->pgfpar[1]);
  parms_SolverSetParam(solver, DTOL,  buf);
}

void fread_param_(char *fname, fprm *prm, int *mprocx, int *mprocy, int
		  *xnmesh, int *ynmesh, int len)
{
  char *mname;

  mname = malloc((len+1)*sizeof(*mname));
  strncpy(mname, fname, len);
  mname[len] = '\0';
  *prm = malloc(sizeof(**prm));
  read_param(mname, *prm);

  *mprocx = (*prm)->mprocx;
  *mprocy = (*prm)->mprocy;
  *xnmesh = (*prm)->xnmesh;
  *ynmesh = (*prm)->ynmesh;

  free(mname);
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
