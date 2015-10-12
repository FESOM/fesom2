/*!
  \file   parms_viewer_impl.h
  \brief  parms_Viewer object.

  \author zzli
  \date   2006-05-05
*/

#ifndef _PARMS_VIEWER_IMPL_H_
#define _PARMS_VIEWER_IMPL_H_

#include <assert.h>
#include "parms_viewer.h"
#include "parms_mem.h"

/*! \struct parms_Viewer_
  \brief parms_Viewer_ structure.
 */
struct parms_Viewer_ {
  int ref;
  BOOL isstd;			//!< is stdout or stderr
  FILE *fp;			//!< file pointer 
  char *fname;			//!< file name 
};

#endif 
