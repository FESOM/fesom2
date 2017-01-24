/*--------------------------------------------------------------------
  parms_PCCreate          : create a preconditioner object.
  parms_PCCreateAbstract  : create an abstract preconditioner object. 
  parms_PCFree            : free the memory for the pc.
  parms_PCILU            : Perform the ilu factorization for the local preconditioner.
  parms_PCSetBsize        : set the block size for ARMS.
  parms_PCSetFill         : set the fill-in parameter for ILUT and
                            ARMS. 
  parms_PCSetILUType      : set the type of ILU.
  parms_PCSetInnerEps     : set the convergence tolerance for the
                            inner GMRES.
  parms_PCSetInnerKSize   : set the restart size for the inner GMRES. 
  parms_PCSetInnerMaxits  : set the maximum iteration counts for the
                            inner GMRES.
  parms_PCSetNlevels      : set the number of levels for ARMS.
  parms_PCSetOP           : set matrix for the preconditioner object.
  parms_PCSetParams       : set parameters for the preconditioner
                            object. 
  parms_PCSetTol          : set the drop tolerance for ILUT and ARMS.
  parms_PCSetTolInd       : set the drop tolerance for finding
                            independent sets in ARMS
  parms_PCSetType         : set the type of preconditioner.
  parms_PCSetup           : set up the preconditioner.
  parms_PCApply           : applying preconditioner.
  parms_PCGetRatio        : get the ratio of the number of nonzeros of 
                            the original matrix to that of the
                            preconditioning matrix.  
  parms_PCView            : dump the preconditioner object. 

  A code fragment for using preconditioning functions:

  parms_PC  pc;
  parms_mat mat;

  // create a preconditioner
  parms_PCCreate(&pc, mat);
  // set the type of the preconditioner.
  parms_PCSetType(pc, PCRAS);
  // set ILU type
  parms_PCSetILUType(pc, PCILUT);
  // set fill-in parameters
  parms_PCSetFill(pc, fill);
  // set drop tolerances
  parms_PCSetTol(pc, tol);
  // set up the preconditioning matrix.
  parms_PCSetup(pc);
  ...
  // free the memory for pc.
  parms_PCFree(&pc);

  $Id: parms_pc.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
  ------------------------------------------------------------------*/
#include "./include/parms_pc_impl.h"
#include "./include/parms_opt_impl.h"

/*------------Protos-------------*/
/*
int parms_ilu0_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
int parms_iluk_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
int parms_ilut_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
int parms_arms_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 

int parms_PCCreate_BJ(parms_PC self);
int parms_PCCreate_Schur(parms_PC self);
int parms_PCCreate_RAS(parms_PC self);
*/
/*----------End Protos------------------*/

/** 
 * Free the memory for the preconditioner object pointed to by self.
 * 
 * @param self A pointer to the memory for the preconditioner object.
 * 
 * @return 0 on success.
 */
int parms_PCFree(parms_PC *self)
{
  (*self)->ref--;
  if ((*self)->ref == 0 ) {
    (*self)->ops->pc_free(self);
    parms_MatFree(&(*self)->A);
    if ((*self)->isperm) {
      PARMS_FREE((*self)->perm);
      PARMS_FREE((*self)->iperm);
    }
    PARMS_FREE((*self)->param);
    PARMS_FREE((*self)->ops);
    PARMS_FREE(*self);
  }

  return 0;
}

/** 
 * Set up the preconditioner (create the preconditioning matrix). 
 * 
 * @param self A preconditioner object.
 * 
 * @return 0 on success.
 */
int parms_PCSetup(parms_PC self)
{
  parms_FactParam param;

  param = self->param;
  
/*--- Define default preconditioner (Block Jacobi/ILU0) ---*/
  if(self->isiluset == false){
    parms_PCSetILUType(self, PCILU0);
  }
  if(self->istypeset == false){
    parms_PCSetType(self, PCBJ);  
  }
/*--- end definition of default precon ---*/
  
  if (self->pctype != PCSCHUR) {
    param->ipar[4] = 0;
    param->ipar[5] = 0;
  }
  
  /* reuse preconditioner (AF)
    if (self->issetup == false) {
    self->issetup = true;
    //self->ops->setup(self);
  }
  else
  {
  // Re-use the precon data struct
     self->issetup = false;
     self->istypeset = false;
     self->ops->pc_free(&self);
     parms_PCSetType(self, self->pctype);
     self->issetup = true;       
     self->ops->setup(self);          
  } */    


  self->issetup = true;
  self->ops->setup(self);
  
  return 0;
}

/** 
 * Create a preconditioner object based on the matrix A.
 * 
 * @param self A preconditioner object.
 * @param A    A matrix object.
 * 
 * @return 0 on success.
 */
int parms_PCCreate(parms_PC *self, parms_Mat A)
{
  parms_PC new_pc;
  int i;

  PARMS_NEW0((new_pc));
  new_pc->ref = 1;
  PARMS_NEW0((new_pc)->ops);
  new_pc->isiluset  = false;
  new_pc->istypeset = false;
  new_pc->issetup   = false;
  new_pc->isopset   = true;
  new_pc->isperm    = false;
  new_pc->A = A;
  A->ref++;
  PARMS_NEW(new_pc->param);
  /* set up default values for  */
  new_pc->param->mc = 1;
  new_pc->param->isalloc = false;
  for (i = 0; i < 7; i++) {
    new_pc->param->lfil[i] = 10;
  }
  new_pc->param->ipar[0] = 5;
  new_pc->param->ipar[1] = 1;
  new_pc->param->ipar[2] = 20;
  new_pc->param->ipar[3] = 0;
  new_pc->param->ipar[4] = 0;
  new_pc->param->ipar[5] = 0;
  for (i = 6; i < 18; i++) {
    new_pc->param->ipar[i] = 0;
  }

  for (i = 0; i < 7; i++) {
    new_pc->param->droptol[i] = 0.001;
  }
  new_pc->param->tolind = 0.05;

  new_pc->param->pgfpar[0] = 0.001;
  new_pc->param->pgfpar[1] = 0.001;

  *self = new_pc;
  return 0;
}

/** 
 * Create an abstract preconditioner object.
 * 
 * @param self A pointer to the preconditioner object.
 * 
 * @return 0 on success.
 */
int parms_PCCreateAbstract(parms_PC *self)
{
  parms_PC new_pc;

  PARMS_NEW(new_pc);
  PARMS_NEW0((new_pc)->ops);
  new_pc->isiluset  = false;
  new_pc->istypeset = false;
  new_pc->issetup   = false;
  new_pc->isopset   = false;
  new_pc->isperm    = false;
  PARMS_NEW(new_pc->param);
  new_pc->param->isalloc = false;

  *self = new_pc;
  return 0;
}

/** 
 * Solve \f$self z = y\f$
 * 
 * @param self A preconditioner object. 
 * @param y    A right-hand-side vector.
 * @param z    The solution vector.
 * 
 * @return 0 on success.
 */
int parms_PCApply(parms_PC self, FLOAT *y, FLOAT *z)
{
  return self->ops->apply(self, y, z);
}

/** 
 * Set the matrix to create the preconditioning matrix.
 * 
 * @param self A preconditioner object.              
 * @param A    The matrix to be used for creating PC.
 * 
 * @return 0 on success.
 */
int parms_PCSetOP(parms_PC self, parms_Mat A)
{

  self->A = A;
  if(!self->isopset){
    A->ref++;
    self->isopset = true;
  }
//  self->issetup = false;
  return 0;
}

/** 
 * Dump preconditioner object self.
 * 
 * @param self A preconditioner object.
 * @param v    A viewer object.
 */
void parms_PCView(parms_PC self, parms_Viewer v)
{
  self->ops->pc_view(self, v);
}

/** 
 * Set preconditioner type.
 *
 *  Currently supported preconditioners:
 *  \f{tabular}{ll}
 *  PCBJ     & block Jacobi \\
 *  PCRAS    & restricted additive Schwarz \\
 *  PCSCHUR  & Schur complement 
 *  \f}
 * 
 * @param self   A preconditioner object.         
 * @param pctype The type of preconditioner.      
 *               - PCBJ    block Jacobi 
 *               - PCRAS   restricted additive Schwarz 
 *               - PCSCHUR  Schur complement 
 * 
 * @return 0 on success.
 */
int parms_PCSetType(parms_PC self, PCTYPE pctype)
{
  parms_Map    is;
  parms_Mat    A;
  BOOL         isserial;

  if (self->istypeset && self->pctype == pctype) {
    return 0;
  }

  A = self->A;
  is = A->is;
  isserial = is->isserial;
  if (self->issetup) {
    self->ops->pc_free(&self);
  }

  if (isserial) {
    self->pctype = PCBJ;
    parms_PCCreate_BJ(self);
    self->istypeset = true;
    self->issetup   = false;
  }
  else {
    self->pctype = pctype;
    if(pctype == PCBJ)
	parms_PCCreate_BJ(self);
    else if(pctype == PCRAS)
	parms_PCCreate_RAS(self);
    else if(pctype == PCSCHUR)
	parms_PCCreate_Schur(self);
    else if(pctype == PCSCHURRAS)
        parms_PCCreate_Schurras(self);
    else
    {
	printf("ERROR: Invalid choice of preconditioner - (Check PCTYPE)! \n");	
	PARMS_ABORT(15);
    }
    self->istypeset = true;
    self->issetup = false;
  }
  return 0;    
}

/** 
 * Set local preconditioner type.
 *
 * Supported ILU preconditioners:
 * 
 *  \f{tabular}{ll}
 *  PCILU0    &  ILU0  \\
 *  PCILUK    &  ILUK  \\
 *  PCILUT    &  ILUT in SPARSKIT \\
 *  PCARMS    &  ARMS implemented by Yousef Saad
 *  \f}
 
 * @param self     A preconditioner object.         
 * @param pcstype  The type of local preconditioner:
 *                 - PCILU0
 *                 - PCILUK
 *                 - PCILUT
 *                 - PCARMS
 *
 * @return 0 on success.
 */
int parms_PCSetILUType(parms_PC self, PCILUTYPE pcilutype)
{

  if (self->isiluset && self->pcilutype == pcilutype) {
    return 0;
  }
  if (self->issetup) {
    self->ops->pc_free(&self);
  }
  self->pcilutype = pcilutype;
  self->A->ilutype = pcilutype;
  self->isiluset = true;
  self->issetup  = false;
  return 0;
}

/* Perform ILU factorization for the local preconditioner*/
int parms_PCILU(parms_PC self, parms_FactParam param, void *mat,
		 parms_Operator *op) 
{
  int type = self->pcilutype;
 
  if(type == PCARMS)
      parms_arms_vcsr(self->A, param, mat, op); 
  else 
      if(*op == NULL){
	  if(type == PCILUT)
	      parms_ilut_vcsr(self->A, param, mat, op); 
	  else if(type == PCILUK)
	      parms_iluk_vcsr(self->A, param, mat, op); 	   
	  else if(type == PCILU0)
	      parms_ilu0_vcsr(self->A, param, mat, op); 
	  else
	  {
	      printf("ERROR: Invalid choice of local preconditioner - (check pcilutype for parms_PCSetILUType(...))\n");
	      PARMS_ABORT(16);
	  }
      }
      else 
	parms_ilu_update(self->A, param, mat, op);
  
  return 0;
}

/** 
 * Set parameters for the preconditioner object.
 *
 * Supported parameters:
 *   - tol    drop tolerance
 *   - fil    fill-in
 *   - nlev   number of levels
 *   - bsize  block size for finding independent sets in ARMS.
 *   - tolind drop tolerance for finding independent sets.
 *   - iksize the restart size for the inner GMRES.
 *   - imax   the number of iterations for the inner GMRES.
 *
 * @param self   A preconditioner object.
 * @param nflags The number of parameters.
 * @param params A pointer to parameters.
 * 
 * @return 0 on success.
 */
int parms_PCSetParams(parms_PC self, int nflags, char **params)
{
  int i, j, k;

  for (i = 0, j = 0; i < nflags; i++) {
    if (!strcmp(params[j], "fill")) {
	for (k = 0; k < 7; k++) { 
	  self->param->lfil[k] = atoi(params[++j]);
	}
    }
    else if (!strcmp(params[j], "tol")) {
      for (k = 0; k < 6; k++) {
	self->param->droptol[k] = atof(params[++j]);
      }
    }
    else if (!strcmp(params[j], "nlev")) {
      self->param->ipar[0] = atoi(params[++j]);
    }
    else if (!strcmp(params[j], "bsize")) {
      self->param->ipar[2] = atoi(params[++j]);
    }
    else if (!strcmp(params[j], "iksize")) {
      self->param->ipar[4] = atoi(params[++j]);
    }
    else if (!strcmp(params[j], "imax")) {
      self->param->ipar[5] = atoi(params[++j]);
    }
    else if (!strcmp(params[j], "tolind")) {
      self->param->tolind = atof(params[++j]);
    }
  }
  return 0;
}

/** 
 * Set fill-in parameter for ILUT and ARMS.
 * 
 * @param self A preconditioner object.                
 * @param fill A int array of size 7.
 *             - fill[0] amount of fill-in kept  in \f$L_{B}\f$. 
 *             - fill[1] amount of fill-in kept  in \f$U_{B}\f$.
 *             - fill[2] amount of fill-in kept  in \f$E L^{-1}\f$.
 *             - fill[3] amount of fill-in kept  in \f$U^{-1}_{B} F\f$.
 *             - fill[4] amount of fill-in kept  in \f$S\f$.
 *             - fill[5] amount of fill-in kept  in \f$L_S\f$.
 *             - fill[6] amount of fill-in kept  in \f$U_S\f$. 
 * 
 * @return 0 on success.
 */
int parms_PCSetFill(parms_PC self, int *fill)
{
  int i;

  for (i = 0; i < 7; i++) {
    self->param->lfil[i] = fill[i];
  }
  return 0;
}

/** 
 * Set the number of levels for ILUK and ARMS.
 * 
 * @param self   A preconditioner object.
 * @param nlevel The number of levels.     
 * 
 * @return 0 on success.
 */
int parms_PCSetNlevels(parms_PC self, int nlev)
{
  self->param->ipar[0] = nlev;
  return 0;
}

/** 
 * Set the type of permutation in ARMS
 * 
 * @param self A preconditioner object.
 * @param type Permutation type.
 *             - 1 non-symmetric permutaion.
 *             - 0 symmetric permutation.
 * 
 * @return 0 on success.
 */
int parms_PCSetPermType(parms_PC self, int type)
{
  self->param->ipar[1] = type;
  return 0;
}

/** 
 * Set the block size for ARMS.
 * 
 * @param self  A preconditioner object.
 * @param bsize The block size for ARMS.
 * 
 * @return 0 on success.
 */
int parms_PCSetBsize(parms_PC self, int bsize)
{
  self->param->ipar[2] = bsize;
  return 0;
}

/** 
 * Set the drop tolerance for ILUT preconditioner.
 * 
 * @param self A preconditioner object.
 * @param tol  A double array of size 7. 
 *             - tol[0]  threshold for dropping  in L_{B}. 
 *             - tol[1]  threshold for dropping  in U_{B}.
 *             - tol[2]  threshold for dropping  in L^{-1} F 
 *             - tol[3]  threshold for dropping  in E U^{-1} 
 *             - tol[4]  threshold for dropping  in Schur complement
 *             - tol[5]  threshold for dropping  in L in last block
 *             - tol[6]  threshold for dropping  in U in last block
 * 
 * @return 0 on success.
 */
int parms_PCSetTol(parms_PC self, double *dt)
{
  PARMS_MEMCPY(self->param->droptol, dt, 7);
  return 0;
}

/** 
 * Set the restart size for the inner GMRES.
 * 
 * @param self A preconditioner object.
 * @param im   The restart size of the inner GMRES.
 * 
 * @return 0 on success.
 */
int parms_PCSetInnerKSize(parms_PC self, int im)
{
  self->param->ipar[4] = im;
  return 0;
}

/** 
 * Set the maximum iteration counts for the inner GMRES.
 * 
 * @param self A preconditioner object.     
 * @param imax The maximum iteration counts.
 * 
 * @return 0 on success.
 */
int parms_PCSetInnerMaxits(parms_PC self, int imaxits)
{
  self->param->ipar[5] = imaxits;
  return 0;
}

/** 
 * Set the convergence tolerance for the inner GMRES.
 * 
 * @param self A preconditioner object.  
 * @param eps  The convergence tolerance.
 * 
 * @return 0 on success.
 */
int parms_PCSetInnerEps(parms_PC self, REAL ieps)
{
  self->param->pgfpar[0] = ieps;
  self->param->pgfpar[1] = ieps;
  return 0;
}

/** 
 * Set the tolerance for finding independent sets.
 * 
 * @param self    A preconditioner object.
 * @param tolind  The drop tolerance for finding independent sets.
 * 
 * @return 0 on success.
 */
int parms_PCSetTolInd(parms_PC self, REAL tolind)
{
  self->param->tolind = tolind;
  return 0;
}

/** 
 * Set permutation and scaling options for interlevel blocks
 * 
 * @param self A preconditioner object.
 * @param meth Options:
 *        -meth[0] nonsummetric permutations of  1: yes. affects rperm
 *                    USED FOR LAST SCHUR COMPLEMENT 
 *        -meth[1] permutations of columns 0:no 1: yes. So far this is
 *                 USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. 
 *                 (so ipar[11] does no matter - enter zero). If 
 *                  ipar[15] is one then ILUTP will be used instead 
 *                  of ILUT. Permutation data stored in: perm2. 
 *         -meth[2] diag. row scaling. 0:no 1:yes. Data: D1
 *         -meth[3] diag. column scaling. 0:no 1:yes. Data: D2
 * @param flag 
 *        -1 interlevel block.
 *        -0 last block.
 * 
 * @return 0 on success.
 */
int parms_PCSetPermScalOptions(parms_PC self, int *meth, int flag) 
{
  if (flag) {
    PARMS_MEMCPY(&self->param->ipar[10], meth, 4);
  }
  else {
    PARMS_MEMCPY(&self->param->ipar[14], meth, 4);
  }

  return 0;
}

/** 
 * Get the ratio of the number of nonzero entries of the
 * preconditioning matrix to that of the original matrix.
 * 
 * @param self  A preconditioner.            
 * @param ratio A pointer to the ratio.      
 * 
 * @return 0 on success.
 */
int parms_PCGetRatio(parms_PC self, double *ratio)
{
  return self->ops->getratio(self, ratio);
}

/** 
 * Return the name of a preconditioner.
 * 
 * @param self A preconditioner.               
 * @param name The name of preconditioner.     
 * 
 * @return 0 on success.
 */
int parms_PCGetName(parms_PC self, char **name)
{
   if(self->pctype == PCBJ)
	*name = "Block Jacobi";
   else if(self->pctype == PCSCHUR)
	*name = "Schur Complement based Preconditioner";
   else if(self->pctype == PCRAS)
	*name = "Restricted Additive Schwarz";   
   else if(self->pctype == PCSCHURRAS)
	*name = "Schur Complement + Restricted Additive Schwarz";   
   else
	*name = "Unknown Preconditioner";
  return 0;
}

/** 
 * Return the name of a local preconditioner.
 * 
 * @param self    A preconditioner.                         
 * @param iluname The name of local ILU preconditioner.     
 * 
 * @return 0 on success.
 */
int parms_PCILUGetName(parms_PC self, char **iluname)
{
   if(self->pcilutype == PCILU0)
	*iluname = "ILU0";
   else if(self->pcilutype == PCILUK)
	*iluname = "ILUK";
   else if(self->pcilutype == PCILUT)
	*iluname = "ILUT";   
   else if(self->pcilutype == PCARMS)
	*iluname = "ARMS";
   else
	*iluname = "Unknown Local Preconditioner";

  return 0;
}
