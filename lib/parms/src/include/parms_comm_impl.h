/*!
  \file   parms_comm_impl.h
  \brief  parms_Comm structure file

  \author zzli
  \date   2006-05-07
*/

#ifndef _PARMS_COMM_IMPL_H_
#define _PARMS_COMM_IMPL_H_

#include "parms_comm.h"
#include "parms_mem.h"

/*! 
  \struct parms_Comm_
 */
struct parms_Comm_ {

  int ref;
  MPI_Comm     comm;		//!< The communicator. 
  COMMTYPE     ctype;		//!< The communication style.
  /*! 
   \param npsend  The number of processors receiving data from the local
    processor.  
   */
  int          npsend; 		
  int          *procs_send;	//!< list of processors to be sent.
  int          *vlist_send;	//!< list of variables to be sent.
  /*! 
   \param ptrvsend An array of pointers to the beginning of the data
   sent to the i-th processor in vlist_send.   
  */
  int          *ptrvsend;	
  FLOAT        *buf_send;	//!< The send buffer.
  MPI_Status   *status_send;	//!< The send status.
  MPI_Request  *req_send;	//!< The send request.
  /*! 
   \param isdt_alloc Indicate if the derived data type for sending
   data is created or not.    
  */
  BOOL         isdt_alloc;
  /*! 
   \param mdata_send The maximum number of the data to be sent on each row.
   */
  int          mdata_send;	
  MPI_Datatype *dtype_send;	//!< The derived send datatype
  /*! 
   \param nprecv The number processors from which data are received
   */
  int          nprecv;
  int          nodv;		//!< The number of external variables.
  /*! 
   \param odvlist The list of external variables stored in processor
   succession. 
   */
  int          *odvlist;	
  /*! 
   \param novp The number of subdomains to which the variable is sent
   */
  int          *novp;
  /*! 
   \param ptrvrecv An array of the beginning of data received from the
   i-th processor.   
  */
  int          *ptrvrecv;
  /*! 
   \param procs_recv The list of processors from which the data are
   received. 
   */
  int          *procs_recv;
  FLOAT        *buf_recv;	 //!< The receive buffer.
  MPI_Status   *status_recv;	//!< The receive status.
  MPI_Request  *req_recv;	//!< The receive request.
};

#endif 
