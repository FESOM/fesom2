/**
 * @file   parms_comm.h
 * @author Zhongze Li
 * @date   Tue Oct 17 10:01:17 2006
 * 
 * @brief  The communication handler used for the matrix-vector
 *         product. 
 * 
 * 
 */

#ifndef _PARMS_COMM_H_
#define _PARMS_COMM_H_

#include "parms_sys.h"
#include "parms_viewer.h"

PARMS_CXX_BEGIN

typedef struct parms_Comm_ *parms_Comm;

/** 
 * Create a parms_Comm object.
 * 
 * @param self A pointer to the parms_Comm object.
 * @param comm MPI communicator
 * 
 * @return 0 on success.
 */
extern int parms_CommCreate(parms_Comm *self, MPI_Comm comm);

/** 
 * Free the memory for the parms_Comm object.
 * 
 * @param self A pointer to the parms_Comm object.
 * 
 * @return 0 on success.
 */
extern int parms_CommFree(parms_Comm *self);

/** 
 * Exchange data for the matrix-vector product.
 *
 * This functions exchanges the interface variables. The offset 
 * indicates the distance between start of the data to be exchanged
 * and the start of the local vector. This is useful for Schur 
 * complement based preconditioner since the data parameter is only
 * the interface part of the local vector rather than the entire local
 * vector. 
 
 * @param self   A parms_Comm object.     
 * @param data   The data to be exchanged.
 * @param offset The distance between the start of the data and the
 * 		 beginning of the local vector.     
 * @return 0 on success.
 */
extern int parms_CommDataBegin(parms_Comm self, void *data, int
			       offset);
/** 
 * Wait for all messages.
 *
 * After calling this function, you may use data in receive buffer
 * safely. 
 *
 * @param self A parms_Comm object.
 * 
 * @return 0 on success.
 */
extern int parms_CommDataEnd(parms_Comm self);

/** 
 * Dump the communication handler comm.
 * 
 * @param self A communication handler. 
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
extern int parms_CommView(parms_Comm self, parms_Viewer v);

/** 
 * Get the total numbjer of variables received.
 * 
 * @param self A communication handler.
 * 
 * @return The number of variables received.
 */
extern int parms_CommGetNumRecv(parms_Comm self);

/** 
 * Get receive buffer.
 * 
 * @param self A communication handler.
 * @param rbuf Receive buffer returned.
 * 
 * @return 0 on success.
 */
extern int parms_CommGetRecvBuf(parms_Comm self, FLOAT **rbuf);

/** 
 * Get the number of variables sent.
 * 
 * @param self A communication handler.
 * 
 * @return The number of vars sent.
 */
extern int parms_CommGetNumSend(parms_Comm self);

PARMS_CXX_END

#endif 
