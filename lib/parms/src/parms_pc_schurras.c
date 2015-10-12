#include "./include/parms_pc_impl.h"
#include "./include/parms_opt_impl.h"
#if defined(__ICC)
#include <mathimf.h>
#else
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#endif 

typedef struct schurras_data {

  parms_Operator op_out,op_in;
  int n, schur_start, nrow, nloc, n_ext, nodv, nsend;
  FLOAT *x_ext, *y_ext;
  BOOL issetup;
  parms_Comm handler;
  parms_Mat S;
  FLOAT *rbuf;
  
} *schurras_data;

static int parms_PC_GetS(parms_PC self, parms_Operator op,parms_Mat *mat);
static int parms_PC_GetExtendSchur(parms_PC self, parms_Operator op,void **mat);
/** Free the memory for struct schur_data.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
static int pc_schurras_free(parms_PC *self)
{
  schurras_data pc_data;
  parms_Operator op;
  int i;

  pc_data = (schurras_data)(*self)->data;
  op = pc_data->op_out;
  parms_OperatorFree(&op);
  op = pc_data->op_in;
  parms_OperatorFree(&op);
  parms_MatFree(&pc_data->S);
 
  if (pc_data->n_ext) {
    PARMS_FREE(pc_data->x_ext);
    PARMS_FREE(pc_data->y_ext);
  }
  PARMS_FREE(pc_data);
  (*self)->param->isalloc = false;    
  return 0;
}

/** Dump the Schur preconditioner.
 *
 *  \param self A preconditioner object.
 *  \param v    A viewer.
 *  \return 0 on success.
 */
static int pc_schurras_view(parms_PC self, parms_Viewer v)
{
  schurras_data pc_data;

  pc_data = (schurras_data)self->data;
  //  parms_OperatorView(pc_data->op_out, v);
  parms_OperatorView(pc_data->op_in, v);
  return 0;
}

/** Set up data for the Schur preconditioner.
 *
 *  Malloc memories for working arrays in the inner GMRES.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
static int pc_schurras_setup(parms_PC self)
{
  schurras_data pc_data;
  parms_Mat A,S;
  parms_PC newpc;
  parms_Viewer v;
  void *diag_mat,*schur_mat,*mat_ext;
  parms_FactParam param;
  parms_Operator op,opras;
  parms_Comm handler;
  parms_Operator newOpt;
  int nrow, i,n_ext,ipar[2];


  A =  self->A;
  pc_data = (schurras_data)self->data;
 
  /* Schur part I */

  /* Perform restricted Gauss Elimination*/
  param = self->param;
   
  param->start = 0;
  param->n = A->is->lsize;  //n = size of the Submatrix
  param->schur_start = A->is->schur_start;  
  param->isalloc = false;
  
  parms_MatGetDiag(A, &diag_mat);  
  op = NULL;
  parms_PCILU(self, param, diag_mat, &op);
  pc_data->op_out = op;
  /* Perform global Schur matrix */
  parms_PC_GetS(self,pc_data->op_out,&S); 
  
  /* rAS-part of the preconditioner */
  
  /* get communicator of global Schur matrix */
  parms_MatGetCommHandler(S,&handler);
  
  /* get Si' = extended local part of global Schur matrix */
  parms_MatGetSubMat(S,&schur_mat);
  parms_MatExtend(S,handler,0,schur_mat,&n_ext,&mat_ext);

  pc_data->handler = handler;
  pc_data->nodv = parms_CommGetNumRecv(handler);
  pc_data->nsend = parms_CommGetNumSend(handler);
  parms_CommGetRecvBuf(handler,&pc_data->rbuf);
    
  /* perform ILU-factorization of Si' */
    
  param->n = n_ext;
  param->schur_start = param->n;
  param->isalloc = false;

  PARMS_NEW(newpc);
  newpc->A = S;
  newpc->pcilutype = self->pcilutype;  
  op = NULL;
  parms_PCILU(newpc,param,mat_ext,&op);
  pc_data->op_in = op;
  PARMS_FREE(newpc);

  /* free the memory for the mat_ext */
  if (self->pcilutype != PCARMS) {
    parms_MatFreeSubMat(S, mat_ext);
  }
 
  /* set pc_data */   
  pc_data->S = S;
  pc_data->n_ext = n_ext;
  pc_data->nrow = nrow = S->is->lsize;
  pc_data->schur_start = parms_OperatorGetSchurPos(pc_data->op_out);
  pc_data->nloc = A->is->lsize;
  pc_data->n = pc_data->nloc + pc_data->nodv;

  if(n_ext){
    PARMS_NEWARRAY(pc_data->x_ext, n_ext);
    PARMS_NEWARRAY(pc_data->y_ext, n_ext);
  }
  pc_data->issetup = true;
  return 0;
}

/** Apply the preconditioner to the vector y.
 *
 *  This preconditioner is actually the lsch_xx preconditioners in old
 *  version of pARMS. 
 *
 *  \param self A preconditioner object.
 *  \param y    A right-hand-side vector.
 *  \param x    Solution vector.
 *  \return 0 on success.
 */
static int pc_schurras_apply(parms_PC self, FLOAT *y, FLOAT *x)
{
  /*--------------------------------------------------------------------
    APPROXIMATE LU-SCHUR LEFT PRECONDITONER
    *--------------------------------------------------------------------
    Schur complement left preconditioner.  This is a preconditioner for 
    the global linear system which is based on solving approximately the
    Schur  complement system. More  precisely, an  approximation to the
    local Schur complement  is obtained (implicitly) in  the form of an
    LU factorization  from the LU  factorization  L_i U_i of  the local
    matrix A_i. Solving with this approximation amounts to simply doing
    the forward and backward solves with the bottom part of L_i and U_i
    only (corresponding to the interface variables).  This is done using
    a special version of the subroutine lusol0 called lusol0_p. Then it
    is possible to  solve  for an  approximate  Schur complement system
    using GMRES  on  the  global  approximate Schur  complement  system
    (which is preconditionied  by  the diagonal  blocks  represented by
    these restricted LU matrices).
    --------------------------------------------------------------------
    Coded by Y. Saad - Aug. 1997 - updated Nov. 14, 1997.
    C version coded by Zhongze Li, Jan. 2001, updated Sept, 17th, 2001
    Revised by Zhongze Li, June 1st, 2006.
    --------------------------------------------------------------------
  */
  schurras_data     pc_data;
  parms_Operator op_in,op_out;
  parms_Map      is;
  parms_Mat      A,S;
  parms_FactParam param;
  MPI_Comm       comm;
  parms_Comm handler;
  int schur_start,nrow, im,i,jj,i1,k,k1,its,j,ii,maxits,incx=1,ierr;
  int if_in_continue,if_out_continue;
  FLOAT *rbuf,*x_ext, *y_ext;
  int nodv, nsend,n,nloc, n_ext;
  REAL *c;
  REAL  eps,eps1,t,ro,tloc;  
  REAL gam;    
  
  static int iters = 0;

  /* retrieve the specific data structure for Schur based preconditioners  */ 
  pc_data = (schurras_data)self->data;

  /* get the matrix A */
  A = self->A;   
  is = A->is;  
  /* get the dimension of the local matrix */
  n = pc_data->n; 

  /* get the beginning location of the Schur complement */
  schur_start = pc_data->schur_start;

  /* get the size of the Schur complement */
  nrow = pc_data->nrow;
 
  /* get tolerance */
  param = self->param;
  eps = param->pgfpar[0];
  eps1 = eps;

  x_ext = pc_data->x_ext;
  y_ext = pc_data->y_ext;
  
  // retrieve the operators 
  op_in = pc_data->op_in;
  op_out = pc_data->op_out;
  
  // get ras-stuff 
  handler = pc_data->handler;
  rbuf = pc_data->rbuf;
  nsend = pc_data->nsend;
  nodv = pc_data->nodv;
  nloc = pc_data->nloc;
  n_ext = pc_data->n_ext;

  S = pc_data->S;
  
  parms_OperatorLsol(op_out, y, x);
 
  ///////////////////////// Vernachlaessigung von Eij /////////////////
  // parms_OperatorInvS(op_out,x,x);
  /////////////////////////////////////////////////////////////////////


  
  ///////////////////////////// Si' mit ILU ///////////////////////
  // matrix-vector product
  if (nsend) {
    parms_CommDataBegin(handler,&x[schur_start], 0);
  }
  PARMS_MEMCPY(y_ext,&x[schur_start],nrow);
  if (nodv) {
    parms_CommDataEnd(handler);
    PARMS_MEMCPY(&y_ext[nrow],rbuf,nodv);//pc_data->n_ext - pc_data->nrow);
  } 
  
  parms_OperatorApply(op_in,y_ext,x_ext);
  
  PARMS_MEMCPY(&x[schur_start],x_ext,nrow);
  // backward sweep 
  
  parms_OperatorAscend(op_out, x, x);
  
  
  return 0;
}

/** Get the ratio of the number of nonzero entries of the
 *  preconditioning matrix to that of the original matrix.
 *
 *  \param self   A preconditioner.
 *  \param ratio  A pointer to the ratio.      
 *  \return 0 on success.
 */
static int pc_schurras_getratio(parms_PC self, double *ratio)
{
  schurras_data pc_data;
  parms_Operator op;
  int nnz_mat, nnz_pc;
  int gnnz_mat, gnnz_pc;

  pc_data = (schurras_data)self->data;
  op = pc_data->op_in;
  parms_OperatorGetNnz(op, &nnz_mat, &nnz_pc);
  MPI_Allreduce(&nnz_mat, &gnnz_mat, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD); 
  MPI_Allreduce(&nnz_pc, &gnnz_pc, 1, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD); 
  *ratio = (double)gnnz_pc/(double)gnnz_mat;  
  return 0;
}

/** Create a Schur-complement based preconditioner.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
int parms_PCCreate_Schurras(parms_PC self)
{
  schurras_data data;

  PARMS_NEW(data);
  data->issetup = false;

  self->data = data;
  self->ops->pc_free  = pc_schurras_free;
  self->ops->pc_view  = pc_schurras_view;
  self->ops->apply    = pc_schurras_apply;
  self->ops->setup    = pc_schurras_setup;
  self->ops->getratio = pc_schurras_getratio;
  
  return 0;
}



static int parms_PC_GetS(parms_PC self, parms_Operator op,parms_Mat *mat)
{

  int i,j,k;
  int n,n_schur,start;
  int npro,pid,nodv,nsend,offset;
  int *cc,*displs,*idx,*idxloc,*counts;
  int cnt1,*nnz,cnt2,newsize;
  int *rowj,*offsetptr,*clp;
  FLOAT *rowa;
  FILE *fptr;
  char fname[15];
  parms_Map is,is_schur;
  parms_Mat S;
  parms_Viewer v;
  void *U_loc,*Offdiag;
  parms_Table newidx;
  parms_vcsr U_vcsr,Offdiag_vcsr;

  MPI_Comm comm;
  parms_Comm handler;

  /* get local Schur Matrix  - it's part of matrix U after restricted Gauss Elimination*/
  parms_OperatorGetU(op,&U_loc);  
  U_vcsr = (parms_vcsr)U_loc;
  n = U_vcsr->n;
  
  /* get entries of Interface-Interface-Submatrix Eij */
  parms_MatGetOffDiag(self->A,&Offdiag);
  Offdiag_vcsr = (parms_vcsr)Offdiag;

  start = parms_OperatorGetSchurPos(op);
  n_schur = n - start;

  /* communicator-data */
  is = self->A->is;
  npro = is->npro;
  pid = is->pid;
  comm = is->comm;
  offset = is->start;

  parms_MatGetCommHandler(self->A,&handler);
  nodv = parms_CommGetNumRecv(handler);
  parms_CommGetOdvlist(handler,&offsetptr);

 
  /* displs is local displacement of global Schur matrix S */
  /* cc[i] contains size of local Schur matrix S_i */
  /* newsize = size of global Schur matrix */
  PARMS_NEWARRAY(cc,npro);
  PARMS_NEWARRAY(displs,npro+1);

  MPI_Allgather(&n_schur,1,MPI_INT,cc,1,MPI_INT,comm);
  displs[0] = 0;
  for(i = 0; i < npro; i++)
    displs[i+1] = displs[i] + cc[i];
  newsize = displs[npro];

  /* idx contains global indices of the Schur matrix components */
  PARMS_NEWARRAY(idx,newsize);
  PARMS_NEWARRAY(idxloc,is->lsize);
  
  parms_MapGetGlobalIndices(is,idxloc);
  for(i = 0; i < n_schur; i++)
      idx[displs[pid]+i] = idxloc[is->iperm[start+i]]-offset;
  MPI_Allgatherv(&idx[displs[pid]],n_schur,MPI_INT,idx,cc,displs,MPI_INT,comm);
  PARMS_FREE(idxloc); 

  /* Create IndexSet-Table of new Ordering */
  parms_TableCreate(&newidx,NULL,newsize);
  for(i = 0; i < newsize; i++)
      parms_TablePut(newidx,idx[i],i);

  /* now idx contains the indices of local rows of the global Schur matrix */
  PARMS_RESIZE(idx, n_schur);
  for(i = 0; i < n_schur; i++)
    idx[i] = displs[pid]+i;
  
  /* create global Schur matrix S */
  parms_MapCreateFromPetsc(&is_schur, n_schur, newsize, comm);
  parms_MatCreate(&S, is_schur);
 
  /*count number of entries of S */
  cnt1 = 0;
  for(i = 0; i < n_schur; i++)
    cnt1 += U_vcsr->nnzrow[start+i] + Offdiag_vcsr->nnzrow[i];   
  
  PARMS_NEWARRAY(nnz,n_schur+1);
  PARMS_NEWARRAY(rowj,cnt1);
  PARMS_NEWARRAY(rowa,cnt1); 

  /* fill global Schur Matrix */
  nnz[0] =0;
  for(i = 0; i < n_schur; i++) {
    cnt1 = U_vcsr->nnzrow[start+i];
    if(cnt1){
      /* entries local Schur matrix */
	PARMS_MEMCPY(&rowa[nnz[i]],U_vcsr->pa[start+i],cnt1);
	for(j = 0; j < cnt1; j++)
	    rowj[nnz[i]+j] = U_vcsr->pj[start+i][j] - start + displs[pid];
    }
    cnt2 = Offdiag_vcsr->nnzrow[i];
    if(cnt2){
      /* entries of Eij */
	PARMS_MEMCPY(&rowa[nnz[i]+cnt1],Offdiag_vcsr->pa[i],cnt2);
	for(j = 0; j < cnt2; j++){
	    clp = parms_TableGet(newidx,offsetptr[Offdiag_vcsr->pj[i][j] - is->lsize]);
	    rowj[nnz[i]+cnt1+j] = *clp;
	}
    }
    nnz[i+1] = nnz[i]+cnt1+cnt2;
  }
  parms_MatSetValues(S,n_schur,idx,nnz,rowj,rowa,INSERT);
  parms_MatSetup(S);
  
  PARMS_FREE(cc);
  PARMS_FREE(displs); 
  PARMS_FREE(rowj);
  PARMS_FREE(rowa);
  PARMS_FREE(nnz);  
  PARMS_FREE(idx);
  parms_TableFree(&newidx);
  parms_MapFree(&is_schur);
  
  /*parms_ViewerCreate(&v,"pschur");
  parms_MatView(S,v);
  parms_ViewerFree(&v);*/
 
  *mat = S;

  return 0;
  
}


