/*--------------------------------------------------------------------
  parms_CommCreate      : create a communication handler
  parms_CommDataBegin   : start data communication
  parms_CommDataEnd     : end data communication
  parms_CommFree        : free the memory for the communication handler
  parms_CommGetNumRecv  : get the number of the data to be received 
  parms_CommGetNumSend  : get the number of the data to be sent
  parms_CommGetRecvBuf  : get the receive buffer
  parms_CommView        : dump the communication handler

  A code fragment for using above functions:

  parms_Comm comm;
  double *data;

  // create a communication handler
  parms_CommCreate(comm, MPI_COMM_WORLD);

  // exchange data
  parms_CommDataBegin(comm, data, 0);
  parms_CommDataEnd(comm);
  
  // free the memory for the comm. handler
  parms_CommFree(comm);

  $Id: parms_comm.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
--------------------------------------------------------------------*/

#include "./include/parms_comm_impl.h"

/** 
 * Create a parms_Comm object.
 * 
 * @param self A pointer to the parms_Comm object.
 * @param comm MPI communicator
 * 
 * @return 0 on success.
 */
int parms_CommCreate(parms_Comm *self, MPI_Comm comm)
{
  parms_Comm comm_handler;
  PARMS_NEW0((comm_handler));
  comm_handler->ref = 1;
  MPI_Comm_dup(comm, &comm_handler->comm);

  /* communication via packed data rather than derived datatype */
  comm_handler->ctype = P2P;

  /* the derived datatype created or not? */
  comm_handler->isdt_alloc = false;

  /* the number of PEs to which this PE sends data */
  comm_handler->npsend = 0;

  /* the ids of the PEs to which this PE sends data */
  comm_handler->procs_send  = NULL;

  /* list of local vars sent to other PEs  */
  comm_handler->vlist_send  = NULL;

  /* the pointer array to the beginning position of the data sent to 
     procs_send[i] in vlist_send */
  comm_handler->ptrvsend    = NULL;

  /* send buffer */
  comm_handler->buf_send    = NULL;

  /* MPI send status  */
  comm_handler->status_send = NULL;

  /* MPI send request */
  comm_handler->req_send    = NULL;

  /* derived datatype */
  comm_handler->dtype_send  = NULL;

  /* the number of PEs from which this PE receives data */
  comm_handler->nprecv      = 0;

  /* the ids of the PEs from which this PE receives data */
  comm_handler->procs_recv  = NULL;

  /* number of data received */
  comm_handler->nodv        = 0;

  /* list of data received  */
  comm_handler->odvlist     = NULL;

  /* the pointer array to the beginning position of the data received from 
     procs_send[i] in odvlist */
  comm_handler->ptrvrecv    = NULL;

  /* receive buffer */
  comm_handler->buf_recv    = NULL; 

  /* MPI receive status  */
  comm_handler->status_recv = NULL;

  /* MPI receive request */
  comm_handler->req_recv    = NULL;

  *self = comm_handler;
  return 0;
}

/** 
 * Free the memory for the parms_Comm object.
 * 
 * @param self A pointer to the parms_Comm object.
 * 
 * @return 0 on success.
 */
int parms_CommFree(parms_Comm *self)
{
  int i;

  //  PARMS_VALIDCOMMPTR(*self);
  /* decrease the reference number */
  (*self)->ref--;
  /* if reference number is zero (the object is not referenced by
     other objects), then remove the object*/
  if ((*self)->ref == 0 ) {
    if ((*self)->npsend) {
      PARMS_FREE((*self)->procs_send);
      PARMS_FREE((*self)->vlist_send);
      PARMS_FREE((*self)->buf_send);
      if ((*self)->isdt_alloc) {
	for (i = 0; i < (*self)->npsend; i++) {
	  MPI_Type_free(&(*self)->dtype_send[i]);
	}
	PARMS_FREE((*self)->dtype_send);
      }
    }
    PARMS_FREE((*self)->ptrvsend);
    PARMS_FREE((*self)->ptrvrecv);
    if ((*self)->nprecv) {
      PARMS_FREE((*self)->procs_recv);
      PARMS_FREE((*self)->odvlist);
      PARMS_FREE((*self)->buf_recv);
    }
    if ((*self)->npsend + (*self)->nprecv) {
      PARMS_FREE((*self)->status_send);
      PARMS_FREE((*self)->req_send);
    }
    MPI_Comm_free(&(*self)->comm);
    PARMS_FREE(*self);
  }
  return 0;
}

/** 
 * Dump the communication handler comm.
 * 
 * @param self A communication handler. 
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
int parms_CommView(parms_Comm self, parms_Viewer v)
{
  FILE *fp;
  int i, j;

  parms_ViewerGetFP(v, &fp);
  fprintf(fp, "=====Dump the content of parms_Comm=====\n");
  fprintf(fp, "npsend = %d\n", self->npsend);
  for (i = 0; i < self->npsend; i++) {
    fprintf(fp, "Sending data to processor%d\n", self->procs_send[i]);
  }
  for (i = 0; i < self->npsend; i++) {
    fprintf(fp, "The data in local indices sent to processor %d\n",
	    self->procs_send[i]); 
    for (j = self->ptrvsend[i]; j < self->ptrvsend[i+1]; j++) {
      fprintf(fp, "vlist_send[%d] = %d\n", j, self->vlist_send[j]);
    }
  }

  fprintf(fp, "nprecv = %d\n", self->nprecv);
  for (i = 0; i < self->nprecv; i++) {
    fprintf(fp, "Receiving data from processor%d\n", self->procs_recv[i]);
  }
  fprintf(fp, "========================================\n");
  parms_ViewerStoreFP(v, fp);
  return 0;
}

/** 
 * Start communication among neighbouring processors.
 *
 * If the self->ctype is set to P2P, then the data are copied to the
 * send buffer, self->buf_send, otherwise, the data are sent directly.
 *
 * The parameter, pos, indicates the offset of the parameter, data, to 
 * the beginning of the local part of the distributed vector and
 * matrix.  
 *
 * @param self A communication handler.                                
 * @param data The data to be sent.                                    
 * @param pos  The distance between the beginning of the parameter
 *             data and the beginning of the local part of the
 *             distirbuted vector.   
 * 
 * @return 0 on success.
 */
int parms_CommDataBegin(parms_Comm self, void *data, int pos)
{
  FLOAT *data_send;
  MPI_Comm comm;
  int i, j, index, start, end, length, tag;

  data_send = (FLOAT *)data;
  comm = self->comm;
  tag = 100;
#if defined(DBL_CMPLX)  
  for (i = 0; i < self->nprecv; i++) {
    start = self->ptrvrecv[i];
    end = self->ptrvrecv[i+1];
    length = end - start;    
    MPI_Irecv(&self->buf_recv[start], length, MPI_CMPLX,
	      self->procs_recv[i], tag, comm, &self->req_recv[i]);
  }
  if (self->ctype == P2P) {
    for (i = 0; i < self->npsend; i++) {
      start  = self->ptrvsend[i];
      end    = self->ptrvsend[i+1];
      length = end - start;
      /* copy data to the send buffer */
      for (j = start; j < end; j++) {
	index = self->vlist_send[j];
	self->buf_send[j] = data_send[index-pos];
      }
      MPI_Isend(&self->buf_send[start], length, MPI_CMPLX,
		self->procs_send[i], tag, comm, &self->req_send[i]); 
    } 
  }
  else if (self->ctype == DERIVED) {
    if (self->isdt_alloc == false) {
      int *blen, *disp_index;

      if (self->npsend) {
	PARMS_NEWARRAY(blen,       self->mdata_send);
	PARMS_NEWARRAY(disp_index, self->mdata_send);
	PARMS_NEWARRAY(self->dtype_send, self->npsend);

      }
      for (i = 0; i < self->npsend; i++) {
	start = self->ptrvsend[i];
	end   = self->ptrvsend[i+1];
	length = end - start;
	for (j = start; j < end; j++) {
	  disp_index[j-start] = self->vlist_send[j]-pos;
	  blen[j-start] = 1;
	}
	MPI_Type_indexed(length, blen, disp_index, MPI_CMPLX, 
			 &self->dtype_send[i]);
	MPI_Type_commit(&self->dtype_send[i]);
      }
      self->isdt_alloc = true;
    }

    for (i = 0; i < self->npsend; i++) {
      MPI_Isend(data_send-pos, 1, self->dtype_send[i],
		self->procs_send[i], tag, comm, &self->req_send[i]); 
    }
  }  
#else
  for (i = 0; i < self->nprecv; i++) {
    start = self->ptrvrecv[i];
    end = self->ptrvrecv[i+1];
    length = end - start;    
    MPI_Irecv(&self->buf_recv[start], length, MPI_DOUBLE,
	      self->procs_recv[i], tag, comm, &self->req_recv[i]);
  }
  if (self->ctype == P2P) {
    for (i = 0; i < self->npsend; i++) {
      start  = self->ptrvsend[i];
      end    = self->ptrvsend[i+1];
      length = end - start;
      /* copy data to the send buffer */
      for (j = start; j < end; j++) {
	index = self->vlist_send[j];
	self->buf_send[j] = data_send[index-pos];
      }
      MPI_Isend(&self->buf_send[start], length, MPI_DOUBLE,
		self->procs_send[i], tag, comm, &self->req_send[i]); 
    } 
  }
  else if (self->ctype == DERIVED) {
    if (self->isdt_alloc == false) {
      int *blen, *disp_index;

      if (self->npsend) {
	PARMS_NEWARRAY(blen,       self->mdata_send);
	PARMS_NEWARRAY(disp_index, self->mdata_send);
	PARMS_NEWARRAY(self->dtype_send, self->npsend);

      }
      for (i = 0; i < self->npsend; i++) {
	start = self->ptrvsend[i];
	end   = self->ptrvsend[i+1];
	length = end - start;
	for (j = start; j < end; j++) {
	  disp_index[j-start] = self->vlist_send[j]-pos;
	  blen[j-start] = 1;
	}
	MPI_Type_indexed(length, blen, disp_index, MPI_DOUBLE, 
			 &self->dtype_send[i]);
	MPI_Type_commit(&self->dtype_send[i]);
      }
      self->isdt_alloc = true;
    }

    for (i = 0; i < self->npsend; i++) {
      MPI_Isend(data_send-pos, 1, self->dtype_send[i],
		self->procs_send[i], tag, comm, &self->req_send[i]); 
    }
  }  
#endif

  return 0;
}

/** 
 * Wait for all communication requests.
 *
 * Data in the receive buffer can be used safely after calling this
 * functions. 
 * 
 * @param self A communication handler.  
 * 
 * @return 0 on success.
 */
int parms_CommDataEnd(parms_Comm self)
{
  int npsend, nprecv;
  MPI_Status  *status_send, *status_recv;
  MPI_Request *req_send, *req_recv;

  npsend = self->npsend;
  nprecv = self->nprecv;
  // status_send = self->status_send;
  // status_recv = self->status_recv;
  req_send = self->req_send;
  req_recv = self->req_recv;

  //  MPI_Waitall(nprecv, req_recv, status_recv);
  // MPI_Waitall(npsend, req_send, status_send);

  MPI_Waitall(nprecv, req_recv, MPI_STATUSES_IGNORE);
  MPI_Waitall(npsend, req_send, MPI_STATUSES_IGNORE);
  
  return 0;
}

/** 
 * Get receive buffer. 
 * 
 * @param self A communication handler.
 * @param rbuf Receive buffer.  
 * 
 * @return 0 on success.
 */
int parms_CommGetRecvBuf(parms_Comm self, FLOAT **rbuf)
{

  *rbuf = self->buf_recv;
  return 0;
}

/** 
 * Get the number of the data received.
 * 
 * @param self A communication handler.
 * 
 * @return The number of data received.
 */
int parms_CommGetNumRecv(parms_Comm self)
{
 
  return self->nodv;
}

/** 
 * Get the number of data to be sent.
 * 
 * @param self A communication handler.
 * 
 * @return The number of data to be sent.
 */
int parms_CommGetNumSend(parms_Comm self)
{

  return self->ptrvsend[self->npsend];
}


/* used by SCHURRAS (AF) */
int parms_CommGetOdvlist(parms_Comm self, int **odvlist)
{
  *odvlist = self->odvlist;
  return 0;
}
