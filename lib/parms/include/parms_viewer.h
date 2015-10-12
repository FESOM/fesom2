/**
 * @file   parms_viewer.h
 * @author Zhongze Li
 * @date   Tue Oct 17 12:01:47 2006
 * 
 * @brief  Functions related to parms_Viewer objects. 
 * 
 */

#ifndef _PARMS_VIEWER_H_
#define _PARMS_VIEWER_H_

#include <stdio.h>
#include <string.h>
#include "parms_sys.h"

PARMS_CXX_BEGIN

typedef struct parms_Viewer_ *parms_Viewer;

/** 
 * Create a parms_Viewer object.
 *
 * If PARMS_COUT and PARMS_CERR are input as fname, they stand for
 * standard output and standard error, respectively. Otherwise, each
 * PE create a file "fnameID.dat". ID stands for ID of PE.
 *
 * @param self  A pointer to the parms_Viewer object.
 * @param fname A file name to store data.
 * 
 * @return 0 on success.
 */
extern int parms_ViewerCreate(parms_Viewer *self, char *fname);

/** 
 * Free the memory of the parms_Viewer object.
 * 
 * @param self A pointer to the parms_Viewer object.
 * 
 * @return 0 on success.
 */
extern int parms_ViewerFree(parms_Viewer *self);

/** 
 * Get a pointer to file pointer. 
 * 
 * @param self A parms_Viewer object.           
 * @param fp   A pointer to the file pointer.   
 * 
 * @return 0 on success.
 */
extern int parms_ViewerGetFP(parms_Viewer self, FILE **fp);

/** 
 * Store fp to the parms_Viewer object.
 * 
 * @param self A parms_Viewer object. 
 * @param fp   A file pointer.     
 * 
 * @return 0 on success.
 */
extern int parms_ViewerStoreFP(parms_Viewer self, FILE *fp);

/** 
 * Retrieve the file name.
 * 
 * @param self  A parms_Viewer object.  
 * @param fname The file name retrieved.
 * 
 * @return 0 on success.
 */
extern int parms_ViewerGetFname(parms_Viewer self, char **fname);

/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_viewercreate_(parms_Viewer *self, char *fname, int *ierr,
			 int len);

extern void parms_viewerfree_(parms_Viewer *self, int *ierr);

extern void parms_viewergetfp_(parms_Viewer *self, FILE **fp, int *ierr);

extern void parms_viewergetfname_(parms_Viewer *self, char **fname, int
			   *ierr);

extern void parms_viewerstorefp_(parms_Viewer *self, FILE *fp, int *ierr);

/*
 *
 * End Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif
