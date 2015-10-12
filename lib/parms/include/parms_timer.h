/**
 * @file   parms_timer.h
 * @author Zhongze Li
 * @date   Tue Oct 17 12:00:25 2006
 * 
 * @brief  Functions related to the parms_Timer object. 
 * 
 * 
 */

#ifndef _PARMS_TIMER_H_
#define _PARMS_TIMER_H_

#include "parms_sys.h"
#include "parms_viewer.h"

PARMS_CXX_BEGIN

typedef struct parms_Timer_ *parms_Timer;

/** 
 * Create a parms_Timer object.
 * 
 * @param self A pointer to the parms_Timer object created.
 */
extern void parms_TimerCreate(parms_Timer *self);
  
/** 
 * Reset the parms_Timer object self.
 * 
 * @param self A parms_Timer object.
 * 
 * @return 0 on success.
 */
extern int parms_TimerReset(parms_Timer self);
  
/** 
 * Reset the elapsed time of self to delay.
 * 
 * @param self  A parms_Timer object. 
 * @param delay Reset the elapsed time to delay seconds.        
 * 
 * @return 0 on success.
 */
extern int parms_TimerResetDelay(parms_Timer self, double delay);
  
/** 
 * Pause the parms_Timer object self.
 * 
 * @param self A parms_Timer object. 
 * 
 * @return 0 on success.
 */
extern int parms_TimerPause(parms_Timer self);

/** 
 * Restart the timer.
 * 
 * @param self A parms_Timer object.
 * 
 * @return 0 on success.
 */
extern int parms_TimerRestart(parms_Timer self);
  
/** 
 * Return The wall-clock time since the last call to
 * parms_TimerReset, parms_TimerResetDelay, parms_TimerRestart.
 * 
 * @param self A parms_Timer object.
 * 
 * @return The elapsed wall-clock time in seconds.
 */
extern double parms_TimerGet(parms_Timer self);

/** 
 * Free the memory for the parms_Timer object
 * 
 * @param self A pointer to the parms_Timer object 
 * 
 * @return 0 on success.
 */
extern int parms_TimerFree(parms_Timer *self);

/** 
 * Dump parms_Timer self via parms_Viewer object v.
 * 
 * @param self A parms_Timer object.    
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
extern int parms_TimerView(parms_Timer self, parms_Viewer v);

/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_timercreate_(parms_Timer *self, int *ierr);

extern void parms_timerfree_(parms_Timer *self, int *ierr);

extern void parms_timerget_(parms_Timer *self, double *t, int *ierr);

extern void parms_timerpause_(parms_Timer *self, int *ierr);

extern void parms_timerreset_(parms_Timer *self, int *ierr);

extern void parms_timerresetdelay_(parms_Timer *self, double *delay, int
			    *ierr);

extern void parms_timerrestart_(parms_Timer *self, int *ierr);

extern void parms_timerview_(parms_Timer *self, parms_Viewer *v, int *ierr);

/*
 *
 * End Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif
