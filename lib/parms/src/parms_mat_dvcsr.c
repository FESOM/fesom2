#include "include/parms_opt_impl.h"
#include "include/parms_comm_impl.h"

/*! \typedef parms_dvcsr
  \brief_dvcsr is a synonym for struct parms_dvcsr.
 */  
/*! \struct parms_dvcsr
  \brief 
 */
typedef struct parms_dvcsr {
  /*! \var diag_mat The diagonal matrix stored in vcsr format.
   */
  parms_vcsr    diag_mat;
 /*!  \var offd_mat The off-diagonal matrix stored in vcsr format.
  */
  parms_vcsr    offd_mat;
  /*! \var mvhandler The parms_Comm object for the matrix-vector product. 
   */
  parms_Comm    mvhandler;	

} *parms_dvcsr;

static int MatVec_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y) 
{
  int lsize, i, j, index, *pj, length;
  parms_dvcsr data;
  parms_Comm  handler;
  parms_vcsr diag_mat, offd_mat;
  parms_Map is;
  FLOAT *pa, *offsetptr, s;
  
  parms_VecPerm(x, self->is);
  parms_VecPerm(y, self->is);

  /* extract diagonal and off-diagonal matrix */
  data     = (parms_dvcsr)self->data;
  handler  = data->mvhandler;
  diag_mat = data->diag_mat;
  offd_mat = data->offd_mat;
  is       = self->is;
  lsize    = is->lsize;

  parms_CommDataBegin(handler, x, 0);

  /* perform local matrix vector product */
  for (i = 0; i < lsize; i++) {
    y[i] = 0.0;
    pj = diag_mat->pj[i];
    pa = diag_mat->pa[i];
    for (j = 0; j < diag_mat->nnzrow[i]; j++) {
      index = pj[j];
      y[i]    += pa[j] * x[index];
    }
  }

  /* off-diagonal part matrix-vector product */
  parms_CommDataEnd(handler);
  offsetptr = handler->buf_recv - lsize;
  for (i = 0; i < is->ninf; i++) {
    pj     = offd_mat->pj[i];
    pa     = offd_mat->pa[i];

    for (j = 0; j < offd_mat->nnzrow[i]; j++) 
      y[i+is->nint]   += pa[j] * offsetptr[pj[j]];

  }
  parms_VecInvPerm(x, self->is);
  parms_VecInvPerm(y, self->is);

  return 0;
}

/** 
 * Setup the matrix
 *
 * This is the most important function in the package. The local
 * variables are divided into two categories: interior variables (not
 * receiving data from other processors) and interface variables
 * (receiving data from other processors). The interior variables are
 * divided into independent variables (not sending data to other
 * processors) and output variables (sending data to other processors)
 * for *asymmetric* pattern. The independent variables are listed first,
 * followd by output variables, and interface variables.  The member
 * variable schur_start is the number of independent
 * variables. i.e., the start position of "Schur" complement in the
 * local linear system.  communication handler for the matrix-vector
 * product is also created.
 * 
 * @param self A matrix object.
 * 
 * @return 0 on success.
 */
static int MatSetup_dvcsr(parms_Mat self)
{
  parms_dvcsr data;
  parms_vcsr  diag_mat = NULL, offd_mat = NULL;
  parms_vcsr aux_data;
  parms_Map  is;
  parms_Comm handler;
  int i, j, *perm, *iperm, mypid, pid, npro, lsize, k1, k2; 
  int *rja=NULL, *tja=NULL, nnz, maxnnz, start, end, nrow;
  int gnodv, cindex, index, gindex;
  int *ginfvars, *ginfv2p, *ginfv2ps, *odperm, *neword, *clp, *info,
    *disp, *pj;     
  int npsendrecv, pos, *procslist=NULL;
  int *ptrvsend, *ptrvrecv=NULL, *cmap=NULL;
  FLOAT *ra=NULL, *tma=NULL, val;
  MPI_Comm comm;

  /* extract information arrays */
  is            = self->is;
  npro          = is->npro;
  mypid         = is->pid;
  self->issetup = true;
  data          = (parms_dvcsr)self->data;
  aux_data      = self->aux_data;
  comm          = is->comm;
  handler       = data->mvhandler;
  lsize         = parms_MapGetLocalSize(is);
  perm          = is->perm;
  iperm         = is->iperm;

  /* free member space */
  if(self->aux_data->space)
      PARMS_FREE(self->aux_data->space);

/* allocate memory for the diagonal matrix */
  if(self->isreset)
  {      
    if(self->resetpattern == SAME_NONZERO_STRUCTURE)
    {
    /* copy data into parallel matrix struct */
      offd_mat = data->offd_mat;
      diag_mat = data->diag_mat;
      /* diagonal part */
      for(i=0; i<diag_mat->n; i++)
      {
        PARMS_MEMCPY(diag_mat->pa[i],  aux_data->pa[i], diag_mat->nnzrow[i]);
      }
      /* off diagonal part */
      for(i=0; i<offd_mat->n; i++)
      {
        index = i+is->nint;
        nrow = diag_mat->nnzrow[index];
        PARMS_MEMCPY(offd_mat->pa[i], &aux_data->pa[index][nrow], offd_mat->nnzrow[i]);
      }
      return 0;
    }
    else
    {
       /* free offd_mat for previous data. No need to free diag_mat 
        * for previous data, since array size (lsize) does not change
        */
        diag_mat = data->diag_mat;
        offd_mat = data->offd_mat;
        if(offd_mat->n){
          PARMS_FREE(offd_mat->nnzrow);
          PARMS_FREE(offd_mat->pa);
          PARMS_FREE(offd_mat->pj);
        }
        PARMS_FREE(offd_mat);
                       
       /* free previous communication handler */
        parms_CommFree(&data->mvhandler);
       /* now create new communication handler */
        parms_Comm new_handler;
        parms_CommCreate(&new_handler, comm);
        handler = data->mvhandler = new_handler;
       /* free processor info and table data */
        if(is->vstable) parms_TableFree(&is->vstable);
        if(is->vsend) PARMS_FREE(is->vsend);    
    }
  }
  else{	
      PARMS_NEW(data->diag_mat);
      diag_mat = data->diag_mat;
      PARMS_NEWARRAY(diag_mat->nnzrow, lsize);
      PARMS_NEWARRAY(diag_mat->pj,     lsize);
      PARMS_NEWARRAY(diag_mat->pa,     lsize);
      diag_mat->n = lsize;
   }

  is->nint = lsize-is->ninf;
  is->schur_start = is->nint;
  k1 = 0;
  k2 = is->nint;

  for (i = 0; i < lsize; i++) {
    if (perm[i] == -1) {
      perm[i] = k1++;
    }
    else {
      perm[i] = k2++;
    }
  }

  /* set up inverse permutation array */    
  for (i = 0; i < lsize; i++) {
    k1 = perm[i];
    iperm[k1] = i;
  }
  is->isperm = true;

  /* permute the rows of the local matrix */
  for (i = 0; i < lsize; i++) {
    if (aux_data->nnzrow[i]) {
      PARMS_RESIZE(aux_data->pj[i], aux_data->nnzrow[i]);
      PARMS_RESIZE(aux_data->pa[i], aux_data->nnzrow[i]);
    }
    k1 = perm[i];
    diag_mat->pj[k1]     = aux_data->pj[i];
    diag_mat->pa[k1]     = aux_data->pa[i];
    diag_mat->nnzrow[k1] = aux_data->nnzrow[i];
  }
  
  PARMS_MEMCPY(aux_data->pj,     diag_mat->pj,     lsize);
  PARMS_MEMCPY(aux_data->pa,     diag_mat->pa,     lsize);
  PARMS_MEMCPY(aux_data->nnzrow, diag_mat->nnzrow, lsize);

  /* permute the column indices of the matrix which are corresponding
     to the interior variables */
  for (i = 0; i < is->nint; i++) {
    rja = aux_data->pj[i];
    nnz = aux_data->nnzrow[i];
    ra  = aux_data->pa[i];
    for (j = 0; j < nnz; j++) {
      rja[j] = perm[rja[j]];	/* in C style */
    }
  }

  /* permute the local column indices of the matrix which are
     corresponding to the interface variables */
  maxnnz = 0;
  for (i = is->nint; i < lsize; i++) {
    rja = aux_data->pj[i];
    nnz = aux_data->nnzrow[i];
    if (maxnnz < nnz) {
      maxnnz =  nnz;
    }
    diag_mat->nnzrow[i] = 0;
    for (j = 0; j < nnz; j++) {
      if (rja[j] >= 0) {
	diag_mat->nnzrow[i]++;
	rja[j] = perm[rja[j]]; /* in C style */
      }
    }
  }

  if (maxnnz) {
    PARMS_NEWARRAY(tja, maxnnz);
    PARMS_NEWARRAY(tma, maxnnz);
  }

  /* allocate memory for the off-diagonal matrix */
  PARMS_NEW(data->offd_mat);
  offd_mat = data->offd_mat;
  offd_mat->n = is->ninf;
  if (is->ninf) {
    PARMS_NEWARRAY(offd_mat->nnzrow, is->ninf);
    PARMS_NEWARRAY(offd_mat->pj,     is->ninf);
    PARMS_NEWARRAY(offd_mat->pa,     is->ninf);
  }
  handler->nodv = parms_TableGetSize(is->table) - lsize;
  if (handler->nodv) {
    PARMS_NEWARRAY(handler->odvlist, handler->nodv);
  }

  /* list column indices of diagonal part first */
  for (i = is->nint; i < lsize; i++) {
    rja    = aux_data->pj[i];
    ra     = aux_data->pa[i];
    nnz    = aux_data->nnzrow[i];
    index  = i-is->nint;
    k1     = 0;
    k2     = diag_mat->nnzrow[i];
    nrow   = diag_mat->nnzrow[i];
    for (j = 0; j < nnz; j++) {
      if (rja[j] >= 0) {
	tja[k1]   = rja[j];
	tma[k1++] = ra[j];
      }
      else {
	tja[k2]   = -rja[j]-1;
	clp = parms_MapGlobalToLocal(is, tja[k2]);
	if (clp != NULL) {
	  handler->odvlist[*clp-lsize] = tja[k2];
	}
	tma[k2++] = ra[j];
      }
    }
    PARMS_MEMCPY(rja, tja, nnz);
    PARMS_MEMCPY(ra,  tma, nnz);
    offd_mat->pj[index]     = &rja[nrow];
    offd_mat->pa[index]     = &ra[nrow];
    offd_mat->nnzrow[index] = nnz-nrow;
  }
  if (maxnnz) {
    PARMS_FREE(tja);
    PARMS_FREE(tma);
  }

  /* sort the column indices of the diagonal matrix in ascending
     order */
  for (i = 0; i < lsize; i++) {
    nnz = diag_mat->nnzrow[i];
    rja = diag_mat->pj[i];
    ra  = diag_mat->pa[i];
    for (j = 0; j < nnz-1; j++) {
      for (k1 = j+1; k1 < nnz; k1++) {
	if (rja[j] > rja[k1]) {
	  index = rja[j];
	  rja[j] = rja[k1];
	  rja[k1] = index;

	  val = ra[j];
	  ra[j] = ra[k1];
	  ra[k1] = val;
	}
      }
    }
  }

  PARMS_NEWARRAY(info, npro);
  PARMS_NEWARRAY(disp, npro+1);
  
  MPI_Allgather(&handler->nodv, 1, MPI_INT, info, 1, MPI_INT, comm);
  disp[0] = 0;
  for (i = 0; i < npro; i++) {
    disp[i+1] = disp[i] + info[i];
  }
  
  /* total number of interface variables */
  gnodv = disp[npro];
  if (gnodv) {
    PARMS_NEWARRAY(ginfvars,  gnodv);
    PARMS_NEWARRAY0(ginfv2p,  gnodv);
    PARMS_NEWARRAY(ginfv2ps,  gnodv);
    MPI_Allgatherv(handler->odvlist, handler->nodv, MPI_INT, ginfvars,
		   info, disp, MPI_INT, comm);

    /* variables sent to other processors */
    for (i = 0; i < npro; i++) {
      if (i != mypid) {
	for (j = disp[i]; j < disp[i+1]; j++) {
	  clp = parms_MapGlobalToLocal(is, ginfvars[j]);
	  if (clp != NULL && *clp < lsize) {
	    ginfv2p[j] = mypid;
	  }
	}
      }
    }

    MPI_Allreduce(ginfv2p, ginfv2ps, gnodv, MPI_INT, MPI_SUM, comm);
    PARMS_FREE(ginfv2p);
    PARMS_NEWARRAY(handler->procs_recv, npro);
    for (i = 0; i < npro; i++) {
      info[i] = 0;
    }
    /* processor IDs from which processor mypid receives data */
    handler->nprecv = 0;	
    for (i = disp[mypid]; i < disp[mypid+1]; i++) {
      pid = ginfv2ps[i];
      if (info[pid] == 0) {
	handler->procs_recv[handler->nprecv++] = pid;
      }
      info[pid]++;
    }

    if (handler->nprecv) {
      PARMS_RESIZE(handler->procs_recv, handler->nprecv);
      
      /* sort neighbouring processors received data from in ascending order */  
      for (i = 0; i < handler->nprecv-1; i++) {
	for (j = i+1; j < handler->nprecv; j++) {
	  if (handler->procs_recv[i] > handler->procs_recv[j]) {
	    pid = handler->procs_recv[i];
	    handler->procs_recv[i] = handler->procs_recv[j];
	    handler->procs_recv[j] = pid;
	  }
	}
      }
    }
    else {
      PARMS_FREE(handler->procs_recv);
    }

    PARMS_NEWARRAY(handler->ptrvrecv, handler->nprecv+1);
    handler->ptrvrecv[0] = 0;
    for (i = 0; i < handler->nprecv; i++) {
      pid = handler->procs_recv[i];
      handler->ptrvrecv[i+1] = handler->ptrvrecv[i] + info[pid];
      info[pid] = i;
    }
    
    if (handler->nodv) {
      PARMS_NEWARRAY(odperm, handler->nodv);
      PARMS_NEWARRAY(handler->buf_recv, handler->nodv);

      /* reorder external variables: the data received from the same 
	 processor are listed consecutively, processor by processor */
      j = disp[mypid];
      for (i = 0; i < handler->nodv; i++) {
	pid = ginfv2ps[i+j];
	index = info[pid];
	odperm[i] = handler->ptrvrecv[index]++;
      }
      for (i = handler->nprecv-1; i > 0; i--) {
	handler->ptrvrecv[i] = handler->ptrvrecv[i-1];
      }
      handler->ptrvrecv[0] = 0;

      /* reorder the external interface variables */
      PARMS_NEWARRAY(neword, handler->nodv);
      for (i = 0; i < handler->nodv; i++) {
	index = odperm[i];
	gindex = handler->odvlist[i];
	parms_TablePut(is->table, gindex, index+lsize);
	neword[index] = gindex;
     }

      PARMS_FREE(handler->odvlist);

      PARMS_FREE(odperm);
      handler->odvlist = neword;
    }
    /* change the global column indices of entries in external
       matrix to the local ones */
    for (i = 0; i < offd_mat->n; i++) {
      nnz = offd_mat->nnzrow[i];
      pj  = offd_mat->pj[i];
      for (j = 0; j < nnz; j++) {
	cindex = pj[j];
	clp  = parms_MapGlobalToLocal(is, cindex);
	pj[j]  = *clp;
      }
    }
    PARMS_FREE(ginfv2ps);
    PARMS_FREE(ginfvars);
  }
  else {
    handler->nprecv = 0;
    PARMS_NEWARRAY(handler->ptrvrecv, 1);
    handler->ptrvrecv[0] = 0;
  }

  /* set up send information */
  MPI_Allgather(&handler->nprecv, 1, MPI_INT, info, 1, MPI_INT, comm);
  disp[0] = 0;
  for (i = 0; i < npro; i++) {
    disp[i+1] = disp[i] + info[i];
  }
  
  if (disp[npro]) {
    PARMS_NEWARRAY(procslist, disp[npro]);
    MPI_Allgatherv(handler->procs_recv, handler->nprecv, MPI_INT, 
		   procslist, info, disp, MPI_INT, comm);
  }
  
  handler->npsend = 0;
  PARMS_NEWARRAY(handler->procs_send, npro);
  for (i = 0; i < npro; i++) {
    if (i != mypid) {
      for (j = disp[i]; j < disp[i+1]; j++) {
	if (procslist[j] == mypid) {
	  handler->procs_send[handler->npsend++] = i;
	  break;
	}
      }
    }
  }
  
  if (disp[npro]) {
    PARMS_FREE(procslist);
  }
  
  if (handler->npsend) {
    PARMS_RESIZE(handler->procs_send, handler->npsend);
  }
  else {
    PARMS_FREE(handler->procs_send);
  }
  
  PARMS_NEWARRAY(handler->ptrvsend, handler->npsend+1);
  handler->ptrvsend[0] = 0;
  ptrvsend = handler->ptrvsend;
  ptrvrecv = handler->ptrvrecv;
  npsendrecv = handler->npsend + handler->nprecv;
  if (npsendrecv) {
    PARMS_NEWARRAY(handler->status_send, npsendrecv);
    PARMS_NEWARRAY(handler->req_send,    npsendrecv);
    handler->status_recv = handler->status_send + handler->npsend;
    handler->req_recv    = handler->req_send + handler->npsend;
    for (i = 0; i < handler->npsend; i++) {
      MPI_Irecv(&disp[i], 1, MPI_INT, handler->procs_send[i], 
		100, comm, &handler->req_send[i]);
    }
    
    for (i = 0; i < handler->nprecv; i++) {
      info[i] = ptrvrecv[i+1] - ptrvrecv[i];
      MPI_Isend(&info[i], 1, MPI_INT, handler->procs_recv[i],
		100, comm, &handler->req_recv[i]);
    }
    
    MPI_Waitall(handler->nprecv, handler->req_recv,
		handler->status_recv); 
    MPI_Waitall(handler->npsend, handler->req_send,
		handler->status_send); 
    
    for (i = 0; i < handler->npsend; i++) {
      ptrvsend[i+1] = ptrvsend[i] + disp[i];
    }
    
    nnz = ptrvsend[handler->npsend];
    if (nnz) {
      PARMS_NEWARRAY(handler->buf_send,   nnz);  
      PARMS_NEWARRAY(handler->vlist_send, nnz);  
    }
    
    handler->mdata_send = 0;    
    if (handler->npsend) {
      for (i = 0; i < handler->npsend; i++) {
	nnz = ptrvsend[i+1] - ptrvsend[i];
	if (handler->mdata_send < nnz) {
	  handler->mdata_send = nnz;
	}
      }
    }

    /* handler->vlist_send stores the variables to be sent to adjacent
       processors in global labeling */
    for (i = 0; i < handler->npsend; i++) {
      pos = ptrvsend[i];
      nnz = ptrvsend[i+1] - ptrvsend[i];
      MPI_Irecv(&handler->vlist_send[pos], nnz, MPI_INT,
		handler->procs_send[i], 200, comm,
		&handler->req_send[i]); 
    }
    
    for (i = 0; i < handler->nprecv; i++) {
      pos = ptrvrecv[i];
      nnz = ptrvrecv[i+1] - ptrvrecv[i];
      MPI_Isend(&handler->odvlist[pos], nnz, MPI_INT,
		handler->procs_recv[i], 200, comm,
		&handler->req_recv[i]);
    }
    
    MPI_Waitall(handler->nprecv, handler->req_recv,
		handler->status_recv);  
    MPI_Waitall(handler->npsend, handler->req_send,
		handler->status_send);

    /* change global labelings  to local ones */
    is->ninf_send = 0;
    parms_TableCreate(&is->vstable, NULL, is->ninf);
    if (is->isperm) {
      for (i = 0; i < handler->npsend; i++) {
	for (j = ptrvsend[i]; j < ptrvsend[i+1]; j++) {
	  cindex = handler->vlist_send[j];
	  clp  = parms_MapGlobalToLocal(is, cindex);
	  index  = perm[*clp];
	  handler->vlist_send[j] = index;
	  clp = parms_TableGet(is->vstable, index);
	  if (clp == NULL) {
	    parms_TablePut(is->vstable, index, is->ninf_send++);
	  }
	}
      }
    }

    /* vsend is an array to store variables sent to other processors */ 
    if (is->ninf_send) {
      PARMS_NEWARRAY(is->vsend, is->ninf_send);
      if (is->nint) {
	PARMS_NEWARRAY(cmap, is->nint);
      }
      for (i = 0; i < is->nint; i++) {
	cmap[i] = -1;
      }

      for (i = 0; i < handler->npsend; i++) {
	for (j = ptrvsend[i]; j < ptrvsend[i+1]; j++) { 
	  index = handler->vlist_send[j];
	  clp = parms_TableGet(is->vstable, index);
	  is->vsend[*clp] = index;

	}
      }
      for (i = 0; i < is->ninf_send; i++) {
	if (is->vsend[i] < is->nint) {
	  cmap[is->vsend[i]] = --is->schur_start;
	}
      }

      /* reorder the local variables in the following order:
       * 1. independent variables which are uncoupled with variables in
       *    adjacent processors.
       * 2. variables sent to adjacent processors but the
       *    corresponding rows are decoupled with variables in
       *    the adjacent processors. (asymmetric pattern)
       * 3. interface variables - coupled with variables in adjacent
       *    processors in the sense of receiving data from
       *    the adjacent processors. 
       */
      if (is->schur_start != is->nint) {
	k1 = 0;
	for (i = 0; i < is->nint; i++) {
	  if (cmap[i] == -1) {
	    cmap[i] = k1++;
	  }
	}
	/* permute column indices of the local matrix */
	for (i = 0; i < lsize; i++) {
	  nnz = diag_mat->nnzrow[i];
	  pj = diag_mat->pj[i];
	  for (j = 0; j < nnz; j++) {
	    if (pj[j] < is->nint) {
	      pj[j] = cmap[pj[j]];
	    }
	  }
	}
	
	/* permute rows which are corresponding to interior variables */
	for (i = 0; i < is->nint; i++) {
	  k1 = cmap[i];
	  diag_mat->nnzrow[k1] = aux_data->nnzrow[i];
	  diag_mat->pj[k1]     = aux_data->pj[i];
	  diag_mat->pa[k1]     = aux_data->pa[i];
	}
	PARMS_MEMCPY(aux_data->nnzrow, diag_mat->nnzrow, is->nint);
	PARMS_MEMCPY(aux_data->pj,     diag_mat->pj,     is->nint);
	PARMS_MEMCPY(aux_data->pa,     diag_mat->pa,     is->nint);
	
	/* update perm and iperm */
	for (i = 0; i < lsize; i++) {
	  if (perm[i] < is->nint) {
	    perm[i] = cmap[perm[i]];
	  }
	}

	for (i = 0; i < lsize; i++) {
	  k1 = perm[i];
	  iperm[k1] = i;
	}

	/* update vlist_send */
	for (i = 0; i < handler->npsend; i++) {
	  start = handler->ptrvsend[i];
	  end   = handler->ptrvsend[i+1];
	  for (j = start; j < end; j++) {
	    if (handler->vlist_send[j] < is->nint) {
	      handler->vlist_send[j] = cmap[handler->vlist_send[j]];
	    }
	  }
	}
	
	/* update vsend */
	for (i = 0; i < is->ninf_send; i++) {
	  if (is->vsend[i] < is->nint) {
	    is->vsend[i] = cmap[is->vsend[i]];
	    parms_TablePut(is->vstable, is->vsend[i], i);
	  }
	}
      }
      if (is->nint) {
	PARMS_FREE(cmap);
      }

    }
    
  }

  PARMS_FREE(info);
  PARMS_FREE(disp);

  return 0;
}  
  

static int MatSetCommType_dvcsr(parms_Mat self, COMMTYPE ctype)
{
  parms_dvcsr data;
  parms_Comm  handler;

  data           = (parms_dvcsr)self->data;
  handler        = data->mvhandler;
  handler->ctype = ctype;
  return 0;
}

static int MatMVPY_dvcsr(parms_Mat self, FLOAT alpha, FLOAT *x,
			 FLOAT beta, FLOAT *y, FLOAT *z) 
{
  int lsize, i, j, index, *pj, length;
  parms_dvcsr data;
  parms_Comm handler;
  parms_vcsr diag_mat, offd_mat;
  parms_Map is;
  FLOAT *pa, *offsetptr, s;

/* no need to permute z here, since it only stores the result */  
  parms_VecPerm(x, self->is);
  parms_VecPerm(y, self->is);

  /* extract diagonal and off-diagonal matrix */
  data     = (parms_dvcsr)self->data;
  handler  = data->mvhandler;
  diag_mat = data->diag_mat;
  offd_mat = data->offd_mat;
  is       = self->is;
  lsize    = is->lsize;

  /* post receive and send actions */
  parms_CommDataBegin(handler, x, 0);  

  /* perform local matrix vector product */
  for (i = 0; i < lsize; i++) {
    s = beta * y[i];
    length = diag_mat->nnzrow[i];
    pj = diag_mat->pj[i];
    pa = diag_mat->pa[i];
    for (j = 0; j < length; j++) {
      index = pj[j];
      s    += alpha * pa[j] * x[index];
    }
    z[i] = s;
  }

  parms_CommDataEnd(handler);

  offsetptr = handler->buf_recv - lsize;
  for (i = 0; i < is->ninf; i++) {
    length = offd_mat->nnzrow[i];
    pj     = offd_mat->pj[i];
    pa     = offd_mat->pa[i];
    s      = z[i+is->nint];
    for (j = 0; j < length; j++) {
      s   += alpha * pa[j] * offsetptr[pj[j]];
    }
    z[i+is->nint] = s;
  }

/* inverse permutation. z now has permutation of y, so has to be inverted as well */
  parms_VecInvPerm(x, self->is);
  parms_VecInvPerm(y, self->is);
  parms_VecInvPerm(z, self->is);
  return 0;
}

static int MatVecSchur_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y, int
			     pos)
{
  parms_dvcsr data;
  parms_Comm  handler;
  parms_Map   is;
  parms_vcsr  offd_mat;
  FLOAT       *offsetptr, *pa, s;
  int i, j, length, *pj, lsize;

  is       = self->is;
  lsize    = is->lsize;
  data     = (parms_dvcsr)self->data;
  handler  = data->mvhandler;
  offd_mat = data->offd_mat;

  if ((x == NULL) || (y == NULL)) {
    return 0;
  }

  parms_CommDataBegin(handler, x, pos);
  for (i = pos; i < is->nint; i++) {
    y[i-pos] = 0.0;
  }
  parms_CommDataEnd(handler);
  offsetptr = handler->buf_recv - lsize;
  for (i = 0; i < is->ninf; i++) {
    length = offd_mat->nnzrow[i];
    pj     = offd_mat->pj[i];
    pa     = offd_mat->pa[i];
    s      = 0.0;
    for (j = 0; j < length; j++) {
      s   += pa[j] * offsetptr[pj[j]];
    }
    y[i+is->nint-pos] = s;
  }
  return 0;
}

static int MatGetDiag_dvcsr(parms_Mat self, void **mat)
{
  parms_dvcsr data;
  parms_vcsr  diag_mat;

  data = (parms_dvcsr)self->data;
  diag_mat = data->diag_mat;
  if (self->ilutype == PCARMS) {
    parms_vcsr diag;
    int i, nnz;

    PARMS_NEW(diag);
    diag->n = diag_mat->n;
    PARMS_NEWARRAY(diag->nnzrow, diag->n);
    PARMS_NEWARRAY(diag->pj,     diag->n);
    PARMS_NEWARRAY(diag->pa,     diag->n);
    for (i = 0; i < diag->n; i++) {
      nnz = diag->nnzrow[i] = diag_mat->nnzrow[i];
      if (nnz) {
	PARMS_NEWARRAY(diag->pj[i], nnz);
	PARMS_NEWARRAY(diag->pa[i], nnz);
	PARMS_MEMCPY(diag->pj[i], diag_mat->pj[i], nnz);
	PARMS_MEMCPY(diag->pa[i], diag_mat->pa[i], nnz);
      }
    }
    *mat = diag;
    return 0;
  }
  *mat = diag_mat;
  return 0;
}

static int MatGetOffDiag_dvcsr(parms_Mat self, void **mat)
{
  parms_dvcsr data;
  parms_vcsr  offd_mat;

  data     = (parms_dvcsr)self->data;
  offd_mat = data->offd_mat;
  *mat     = offd_mat;
  return 0;

}

static int MatGetSubMat_dvcsr(parms_Mat self, void **mat)
{
  *mat = self->aux_data;
  return 0;
}

static int MatExtend_dvcsr(parms_Mat self, parms_Comm handler, int
			   start, void *mat, int *nn, void
			   **ext_mat)  
{
  parms_vcsr lmat, new_mat;
  parms_Map  is;
  MPI_Request *sreq=NULL, *rreq=NULL;
  MPI_Status *sstatus=NULL, *rstatus=NULL;
  MPI_Comm comm;
  FLOAT **a_send=NULL, *pa_rbuf=NULL, *pa_sbuf=NULL, *rowa=NULL;
  int **ja_send=NULL;
  int *ia_send=NULL, *pia_send=NULL, *pia_recv=NULL, *pj_rbuf=NULL, *pj_sbuf=NULL;
  int *pindex, *rowj, *rbuf=NULL, *sbuf=NULL;
  int lsize, jcol, i, j, index, nnz, nodv, numsend, numrecv;
  int n, pos, lindex, tag, length, k, m;

  lmat = (parms_vcsr)mat;
  PARMS_NEW(new_mat);
  is = self->is;
  lsize = is->lsize;

  MPI_Comm_dup(handler->comm, &comm);

  if (is->ninf_send) {
    PARMS_NEWARRAY(a_send,  is->ninf_send);
    PARMS_NEWARRAY(ja_send, is->ninf_send);
    PARMS_NEWARRAY(ia_send, is->ninf_send);

    /* copy the matrix corresponding to the interface variables to
       working arrays. the column indices are stored in global numbering */ 
    for (i = 0; i < is->ninf_send; i++) {
      index = is->vsend[i] - start;
      nnz = ia_send[i] = lmat->nnzrow[index];
      PARMS_NEWARRAY(a_send[i],  nnz);
      PARMS_NEWARRAY(ja_send[i], nnz);

      nnz = lmat->nnzrow[index];
      for (j = 0; j < nnz; j++) {
	a_send[i][j] = lmat->pa[index][j];
	jcol = lmat->pj[index][j];
	if (jcol < lsize) {
	  if (is->isperm) {
	    jcol = is->iperm[jcol];
	    ja_send[i][j] = is->lvars[jcol];
	  }
	  else {
	    ja_send[i][j] = is->lvars[jcol];
	  }
	}
	else {
	  jcol = handler->odvlist[jcol-lsize];
	  ja_send[i][j] = jcol;
	}
      }
    }
  }


  numsend = handler->ptrvsend[handler->npsend];
  numrecv = handler->ptrvrecv[handler->nprecv];

  nodv = handler->nodv;
  n = lmat->n + nodv;

  new_mat->n = n;
  PARMS_NEWARRAY(new_mat->nnzrow, n);
  PARMS_NEWARRAY(new_mat->pj,     n);
  PARMS_NEWARRAY(new_mat->pa,     n);

  /* copy the local matrix to new_mat */
  for (i = 0; i < lmat->n; i++) {
    nnz = new_mat->nnzrow[i] = lmat->nnzrow[i];
    if (nnz) {
      PARMS_NEWARRAY(new_mat->pj[i], nnz);
      PARMS_NEWARRAY(new_mat->pa[i], nnz);
      PARMS_MEMCPY(new_mat->pj[i], lmat->pj[i], nnz);
      PARMS_MEMCPY(new_mat->pa[i], lmat->pa[i], nnz);
    }
  }

  if (handler->npsend+handler->nprecv) {
    PARMS_NEWARRAY(sstatus, handler->npsend+handler->nprecv);
    PARMS_NEWARRAY(sreq,    handler->npsend+handler->nprecv);
    rstatus = sstatus + handler->npsend;
    rreq    = sreq + handler->npsend;
  }

  tag = 800;

  /* exchange external matrix */
  if (numrecv) {
    PARMS_NEWARRAY(rbuf, numrecv);
  }

  /* receive the number of entries for each row in the extended
     matrix */
  for (i = 0; i < handler->nprecv; i++) {
    pos = handler->ptrvrecv[i];
    length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
    MPI_Irecv(&rbuf[pos], length, MPI_INT, handler->procs_recv[i],
	      tag, comm, &rreq[i]); 
  }    

  if (numsend) {
    PARMS_NEWARRAY(sbuf, numsend);
  }

  for (i = 0; i < handler->npsend; i++) {
    pos = handler->ptrvsend[i];
    length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
    for (j = 0; j < length; j++) {
      lindex = handler->vlist_send[pos+j];
      pindex = parms_TableGet(is->vstable, lindex);
      lindex = *pindex;
      sbuf[pos+j] = ia_send[lindex];
    }
    MPI_Isend(&sbuf[pos], length, MPI_INT, handler->procs_send[i],
	      tag, comm, &sreq[i]);
  }
  MPI_Waitall(handler->nprecv, rreq, rstatus);
  MPI_Waitall(handler->npsend, sreq, sstatus);

  for (i = 0; i < handler->nprecv; i++) {
    pos = handler->ptrvrecv[i];
    length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];    
    for (j = 0; j < length; j++) {
      new_mat->nnzrow[lmat->n+pos+j] = rbuf[pos+j];
    }
  }
  
  j = 0;
  for (i = 0; i < nodv; i++) {
    length = new_mat->nnzrow[i+lsize];
    j += length;
    if (length) {
      PARMS_NEWARRAY(new_mat->pj[i+lmat->n], length);
      PARMS_NEWARRAY(new_mat->pa[i+lmat->n], length);
    }
  }
  if (j) {
    PARMS_NEWARRAY(pj_rbuf, j);
    PARMS_NEWARRAY(pa_rbuf, j);
  }

  PARMS_NEWARRAY(pia_send, handler->npsend+1);
  PARMS_NEWARRAY(pia_recv, handler->nprecv+1);
  pia_send[0] = 0;
  pia_recv[0] = 0;

  /* receive column indices of each row */
  for (i = 0; i < handler->nprecv; i++) {
    pos = handler->ptrvrecv[i];
    length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
    k = 0;
    for (j = 0; j < length; j++) {
      m = new_mat->nnzrow[lmat->n+pos+j];
      k += m;
    }
    pia_recv[i+1] = pia_recv[i] + k;
    MPI_Irecv(&pj_rbuf[pia_recv[i]],  k, MPI_INT,
	      handler->procs_recv[i], tag, comm, &rreq[i]);    
  }

  k = 0;
  for (i = 0; i < handler->npsend; i++) {
    pos = handler->ptrvsend[i];
    length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
    for (j = 0; j < length; j++) {
      lindex = handler->vlist_send[pos+j];
      pindex = parms_TableGet(is->vstable, lindex);
      lindex = *pindex;
      m = ia_send[lindex];
      k += m;
    }
  }

  if (k) {
    PARMS_NEWARRAY(pj_sbuf, k);
    PARMS_NEWARRAY(pa_sbuf, k);
  }

  for (i = 0; i < handler->npsend; i++) {
    pos = handler->ptrvsend[i];
    length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
    k = 0;
    for (j = 0; j < length; j++) {
      lindex = handler->vlist_send[pos+j];
      pindex = parms_TableGet(is->vstable, lindex);
      lindex = *pindex;
      m = ia_send[lindex];
      PARMS_MEMCPY(&pj_sbuf[pia_send[i]+k], ja_send[lindex], m);
      k += m;
    }
    pia_send[i+1] = pia_send[i] + k;
    MPI_Isend(&pj_sbuf[pia_send[i]], k, MPI_INT, handler->procs_send[i], 
	      tag, comm, &sreq[i]);
  }
  MPI_Waitall(handler->nprecv, rreq, rstatus);
  MPI_Waitall(handler->npsend, sreq, sstatus);

  for (i = 0; i < handler->nprecv; i++) {
    pos = handler->ptrvrecv[i];
    length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
    k = 0;
    for (j = 0; j < length; j++) {
      m = new_mat->nnzrow[lmat->n+pos+j];
      PARMS_MEMCPY(new_mat->pj[lmat->n+pos+j], &pj_rbuf[pia_recv[i]+k],
		   m); 
      k += m;
    }
  }

  /* receive values of each row */

#if defined(DBL_CMPLX)  
  for (i = 0; i < handler->nprecv; i++) {
    k = pia_recv[i+1] - pia_recv[i];
    MPI_Irecv(&pa_rbuf[pia_recv[i]], k, MPI_CMPLX,
	      handler->procs_recv[i], tag, comm, &rreq[i]);  
  }
  
  for (i = 0; i < handler->npsend; i++) {
    pos = handler->ptrvsend[i];
    length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
    k = 0;
    for (j = 0; j < length; j++) {
      lindex = handler->vlist_send[pos+j];
      pindex = parms_TableGet(is->vstable, lindex);
      lindex = *pindex;
      m = ia_send[lindex];
      PARMS_MEMCPY(&pa_sbuf[pia_send[i]+k], a_send[lindex], m);
      k += m;
    }
    MPI_Isend(&pa_sbuf[pia_send[i]], k, MPI_CMPLX,
	      handler->procs_send[i], tag, comm, &sreq[i]);
  }  
#else
  for (i = 0; i < handler->nprecv; i++) {
    k = pia_recv[i+1] - pia_recv[i];
    MPI_Irecv(&pa_rbuf[pia_recv[i]], k, MPI_DOUBLE,
	      handler->procs_recv[i], tag, comm, &rreq[i]);  
  }

  for (i = 0; i < handler->npsend; i++) {
    pos = handler->ptrvsend[i];
    length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
    k = 0;
    for (j = 0; j < length; j++) {
      lindex = handler->vlist_send[pos+j];
      pindex = parms_TableGet(is->vstable, lindex);
      lindex = *pindex;
      m = ia_send[lindex];
      PARMS_MEMCPY(&pa_sbuf[pia_send[i]+k], a_send[lindex], m);
      k += m;
    }
    MPI_Isend(&pa_sbuf[pia_send[i]], k, MPI_DOUBLE,
	      handler->procs_send[i], tag, comm, &sreq[i]);
  }
#endif

  MPI_Waitall(handler->nprecv, rreq, rstatus);
  MPI_Waitall(handler->npsend, sreq, sstatus);

  MPI_Comm_free(&comm);
  for (i = 0; i < handler->nprecv; i++) {
    pos = handler->ptrvrecv[i];
    length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
    k = 0;
    for (j = 0; j < length; j++) {
      m = new_mat->nnzrow[lmat->n+pos+j];
      PARMS_MEMCPY(new_mat->pa[lmat->n+pos+j], &pa_rbuf[pia_recv[i]+k],
		   m); 
      k += m;
    }
  }

  /* free temporary arrays */
  if (handler->npsend) {
    PARMS_FREE(pj_sbuf);
    PARMS_FREE(pa_sbuf);
  }

  if (handler->nprecv) {
    PARMS_FREE(pj_rbuf);
    PARMS_FREE(pa_rbuf);
  }

  if (numrecv) {
    PARMS_FREE(rbuf);
  }

  if (numsend) {
    PARMS_FREE(sbuf);
  }


  PARMS_FREE(pia_recv);
  PARMS_FREE(pia_send);

  if (handler->npsend+handler->nprecv) {
    PARMS_FREE(sreq);
    PARMS_FREE(sstatus);
  }

  if (is->ninf_send) {
    for (i = 0; i < is->ninf_send; i++) {
      if (ia_send[i]) {
	PARMS_FREE(ja_send[i]);
	PARMS_FREE(a_send[i]);
      }
    }

    PARMS_FREE(ia_send);
    PARMS_FREE(ja_send);
    PARMS_FREE(a_send);
  }

  /* change the global column indices into local ones */     
  for (i = lmat->n; i < n; i++) {
    nnz = new_mat->nnzrow[i];
    rowj = new_mat->pj[i];
    rowa = new_mat->pa[i];
    k = 0;
    for (j = 0; j < nnz; j++) {
      jcol = rowj[j];
      pindex = parms_TableGet(is->table, jcol);
      if (pindex) {
	if (*pindex < lsize) {
	  rowj[j] = is->perm[*pindex] - start;
	}
	else if (*pindex >= lsize) {
	  rowj[j] = *pindex - start;
	}
	if (k < j) {
	  rowj[k] = rowj[j];
	  rowa[k] = rowa[j];
	}
	k++;
      }
    }
    new_mat->nnzrow[i] = k;
    if (k) {
      PARMS_RESIZE(new_mat->pj[i], k);
      PARMS_RESIZE(new_mat->pa[i], k);
    }
  }
  *ext_mat = new_mat;
  *nn = n;
  return 0;
}

static int MatFreeSubMat_dvcsr(parms_Mat self, void *mat)
{
  parms_vcsr vmat;
  int i, nnz, isalloc;

  isalloc = self->isalloc;
  
  if(isalloc){
   vmat = (parms_vcsr)mat;
   for (i = 0; i < vmat->n; i++) {
     nnz = vmat->nnzrow[i];
     if (nnz) {
       PARMS_FREE(vmat->pj[i]);
       PARMS_FREE(vmat->pa[i]);
     }
   }

   if (vmat->n) {
     PARMS_FREE(vmat->nnzrow);
     PARMS_FREE(vmat->pj);
     PARMS_FREE(vmat->pa);
   }

   PARMS_FREE(vmat);

   mat = NULL;
  }

  return 0;
}

static int MatGetHandler_dvcsr(parms_Mat self, parms_Comm *handler)
{
  parms_dvcsr data;

  data = (parms_dvcsr)self->data;
  *handler = data->mvhandler;
  return 0;
}

static struct parms_Mat_ops parms_matops_dvcsr = {
  MatVec_dvcsr,
  MatSetup_dvcsr,
  MatSetCommType_dvcsr,
  MatMVPY_dvcsr,
  MatGetDiag_dvcsr,
  MatGetOffDiag_dvcsr,
  MatGetSubMat_dvcsr,
  MatExtend_dvcsr,
  MatVecSchur_dvcsr,
  MatFreeSubMat_dvcsr,
  MatGetHandler_dvcsr
};

int parms_MatFree_dvcsr(parms_Mat *self)
{
  parms_dvcsr   data;
  parms_vcsr    diag_mat, offd_mat;

  data = (parms_dvcsr)(*self)->data;

  parms_CommFree(&data->mvhandler);
  diag_mat = data->diag_mat;
  PARMS_FREE(diag_mat->nnzrow);
  PARMS_FREE(diag_mat->pa);
  PARMS_FREE(diag_mat->pj);
  PARMS_FREE(diag_mat);

  offd_mat = data->offd_mat;
  if (offd_mat->n) {
    PARMS_FREE(offd_mat->nnzrow);
    PARMS_FREE(offd_mat->pa);
    PARMS_FREE(offd_mat->pj);
  }
  PARMS_FREE(offd_mat);

  PARMS_FREE(data);  
  return 0;
}

int parms_MatView_dvcsr(parms_Mat self, parms_Viewer v)
{
  int i, j, lsize, pid, npro, length;
  int *rja, *lvlist, *iperm, index;
  FLOAT *ra;
  parms_dvcsr data;
  parms_vcsr diag_mat, offd_mat;
  parms_Map is;
  FILE *fp;

  is      = self->is;
  lvlist  = is->lvars;
  iperm   = is->iperm;
  lsize   = is->lsize;
  pid     = is->pid;
  npro    = is->npro;
  data    = (parms_dvcsr)self->data;

  parms_ViewerGetFP(v, &fp);
  
  if (pid == 0) {
    fprintf(fp, "There are %d processors available\n", npro);
    fprintf(fp, "The format of output local equations is:\n");
    fprintf(fp, "(pid,local_row_index)=>(pid,global_row_index)\n");
    fprintf(fp, "(pid,local_row_index, local_column_index) = value\n"); 
  }
  fprintf(fp, "\n");

  for (i = 0; i < lsize; i++) {
    index = iperm[i];
    fprintf(fp, "(%d,%d)=>(%d,%d)\n", pid, i, pid, lvlist[index]);
  }
  fprintf(fp, "\n");

  diag_mat = data->diag_mat;
  offd_mat = data->offd_mat;
  fprintf(fp, "Local diagonal matrix on processor %d is\n", pid);

#if defined(DBL_CMPLX)  
  for (i = 0; i < lsize; i++) {
    length = diag_mat->nnzrow[i];
    rja    = diag_mat->pj[i];
    ra     = diag_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "(%d,%d,%d) = (%f, %f)  ", pid, i, rja[j], creal(ra[j]), cimag(ra[j]));
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "Local off-diagonal matrix on processor %d is \n", pid); 
  for (i = 0; i < is->ninf; i++) {
    length = offd_mat->nnzrow[i];
    rja    = offd_mat->pj[i];
    ra     = offd_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "(%d,%d,%d) = (%f, %f)  ", pid, i, rja[j], creal(ra[j]), cimag(ra[j]));
    }
    fprintf(fp, "\n");
  }
#else
  for (i = 0; i < lsize; i++) {
    length = diag_mat->nnzrow[i];
    rja    = diag_mat->pj[i];
    ra     = diag_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "(%d,%d,%d) = %f  ", pid, i, rja[j], ra[j]);
    }
    fprintf(fp, "\n");
  }

  fprintf(fp, "Local off-diagonal matrix on processor %d is \n", pid); 
  for (i = 0; i < is->ninf; i++) {
    length = offd_mat->nnzrow[i];
    rja    = offd_mat->pj[i];
    ra     = offd_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "(%d,%d,%d) = %f  ", pid, i, rja[j], ra[j]);
    }
    fprintf(fp, "\n");
  }
#endif
  fprintf(fp, "\n");
  parms_ViewerStoreFP(v, fp);
  parms_CommView(data->mvhandler, v);
  return 0;
}


int parms_MatViewCOO_dvcsr(parms_Mat self, parms_Viewer v)
{
  int i, j, lsize, pid, npro, length;
  int *rja;

  FLOAT *ra;
  parms_dvcsr data;
  parms_vcsr diag_mat, offd_mat;
  parms_Map is;
  FILE *fp;

  is      = self->is;
  lsize   = is->lsize;
  pid     = is->pid;
  npro    = is->npro;
  data    = (parms_dvcsr)self->data;

  parms_ViewerGetFP(v, &fp);
  
  if (pid == 0) {
    fprintf(fp, "There are %d processors available\n", npro);
    fprintf(fp, "The format of output local equations is:\n");
    fprintf(fp, "local_row_index  local_column_index   value\n"); 
  }
  
  fprintf(fp, "\n");

  diag_mat = data->diag_mat;
  offd_mat = data->offd_mat;

#if defined(DBL_CMPLX)
  for (i = 0; i < lsize; i++) {
    length = diag_mat->nnzrow[i];
    rja    = diag_mat->pj[i];
    ra     = diag_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "%d  %d  %f  %f  \n", i, rja[j], creal(ra[j]), cimag(ra[j]));
    }
  }

  for (i = 0; i < is->ninf; i++) {
    length = offd_mat->nnzrow[i];
    rja    = offd_mat->pj[i];
    ra     = offd_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "%d  %d  %f  %f  \n", i+lsize-is->ninf, rja[j], creal(ra[j]), cimag(ra[j]));
    }
  }
#else
  for (i = 0; i < lsize; i++) {
    length = diag_mat->nnzrow[i];
    rja    = diag_mat->pj[i];
    ra     = diag_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "%d  %d  %f  \n", i, rja[j], ra[j]);
    }
  }

  for (i = 0; i < is->ninf; i++) {
    length = offd_mat->nnzrow[i];
    rja    = offd_mat->pj[i];
    ra     = offd_mat->pa[i];
    for (j = 0; j < length; j++) {
      fprintf(fp, "%d  %d  %f  \n", i+lsize-is->ninf, rja[j], ra[j]);
    }
  }
#endif
  parms_ViewerStoreFP(v, fp);
  return 0;
}


int parms_MatCreate_dvcsr(parms_Mat self)
{
  parms_dvcsr  data;
  parms_Comm   mvhandler;
  MPI_Comm     comm;

  PARMS_MEMCPY(self->ops, &parms_matops_dvcsr, 1);
  PARMS_NEW(data);
  comm = self->is->comm;
  parms_CommCreate(&mvhandler, comm);
  data->mvhandler = mvhandler;
  
  self->data = data;
  return 0;
}
