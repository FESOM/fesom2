/*!
  \file   parms_timer_impl.h
  \brief  parms_Timer structure

  \author zzli
  \date   2006-05-05
*/

#ifndef _PARMS_TIMER_IMPL_H_
#define _PARMS_TIMER_IMPL_H_

#include "parms_mem.h"
#include "parms_timer.h"

/*! \struct parms_Timer_
  \brief parms_Timer_ structure.
 */
struct parms_Timer_ {

  int ref;
  /*! \var initial_time
    \brief initial time
  */
  double initial_time;
  /*! \var elapsed_time
    \brief elapsed time so far.
   */
  double elapsed_time;
};

#endif 
