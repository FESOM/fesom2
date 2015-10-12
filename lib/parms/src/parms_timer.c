/*--------------------------------------------------------------------
  dwalltime              : return the ellapsed time in seconds since
                           the Epoch. 
  parms_TimerCreate      : create a parms_Timer object.
  parms_TimerFree        : free the memory for the parms_Timer object.
  parms_TimerGet         : get the elapsed time since the last call
                           to parms_TimerReset, parms_TimerResetDelay,
                           or parms_TimerRestart.
  parms_TimerPause       : pause the timer.
  parms_TimerRestart     : restart the timer.
  parms_TimerReset       : reset the timer to 0.
  parms_TimerResetDelay  : suspend the timer. 
  parms_TimerView        : dump the parms_Timer object.

  A code fragment for using timer functions:

  parms_Timer t;

  // create a timer.
  parms_TimerCreate(&t);
  parms_TimerReset(t);
  ... fragment code 1
  // time spent on the fragment code 1.
  parms_TimerGet(t);  
  parms_TimerPause(t);
  ...
  ...
  parms_TimerRestart(t);
  ... framgment code 2
  // total time spent on fragment code1 and code2.
  parms_TimerGet(t) 
  // free the memory for t
  parms_TimerFree(&t);
     
  $Id: parms_timer.c,v 1.1.1.1 2006-11-27 22:28:01 zzli Exp $
  ------------------------------------------------------------------*/
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <stdio.h>
#include "./include/parms_timer_impl.h"
#include "parms_viewer.h"

/** 
 * Return the elapsed time since the Epoch.
 * 
 * @return The wall-clock time in seconds.
 */
static double dwalltime()
{
#ifdef USE_MPI
  int flag;
  MPI_Initialized(&flag);
  if (flag) {
    return MPI_Wtime();
  }
  else {
    struct timeval tval;
    double t;

    gettimeofday(&tval, NULL);
    t = (double)(tval.tv_sec + tval.tv_usec  / 1.0e6);
    return t;
  }
#else
  if (true) {
    struct timeval tval;
    double t;
  
    gettimeofday(&tval, NULL);
    t = (double)tval.tv_sec + (double)tval.tv_usec/1.0e6;
    return t;
  }
#endif 
}

/** 
 * Reset the timer to 0.
 *
 * Set the elapsed time to zero and set the initial time of the
 * timer. 
 * 
 * @param self A parms_Timer object.
 * 
 * @return 0 on success.
 */
int parms_TimerReset(parms_Timer self)
{
  self->elapsed_time = 0.0;
  self->initial_time = dwalltime();
  return 0;
}

/** 
 * Create a parms_Timer object.
 * 
 * @param self A pointer to the parms_Timer object created.
 */
void parms_TimerCreate(parms_Timer *self)
{
  parms_Timer new_timer;

  PARMS_NEW0((new_timer));
  new_timer->ref = 1;
  parms_TimerReset(new_timer);
  *self = new_timer;
}

/** 
 * Reset the elapsed time of self to delay.
 * 
 * @param self  A parms_Timer object. 
 * @param delay Reset the elapsed time to delay seconds.        
 * 
 * @return 0 on success.
 */
int parms_TimerResetDelay(parms_Timer self, double delay) 
{
  parms_TimerReset(self);
  self->elapsed_time = delay;
  return 0;
}

/** 
 * Pause the parms_Timer object self.
 * 
 * @param self A parms_Timer object. 
 * 
 * @return 0 on success.
 */
int parms_TimerPause(parms_Timer self)
{
  self->elapsed_time = parms_TimerGet(self);
  return 0;
}

/** 
 * Restart the timer.
 * 
 * @param self A parms_Timer object.
 * 
 * @return 0 on success.
 */
int parms_TimerRestart(parms_Timer self)
{
  self->initial_time = dwalltime();
  return 0;
}

/** 
 * Return The wall-clock time since the last call to
 * parms_TimerReset, parms_TimerResetDelay, parms_TimerRestart.
 * 
 * @param self A parms_Timer object.
 * 
 * @return The elapsed wall-clock time in seconds.
 */
double parms_TimerGet(parms_Timer self)
{
  double current, time;

  current = dwalltime();
  time = current - self->initial_time;
  return (time + self->elapsed_time);
}

/** 
 * Free the memory for the parms_Timer object
 * 
 * @param self A pointer to the parms_Timer object 
 * 
 * @return 0 on success.
 */
int parms_TimerFree(parms_Timer *self)
{

  PARMS_FREE(*self);
  return 0;
}

/** 
 * Dump parms_Timer self via parms_Viewer object v.
 * 
 * @param self A parms_Timer object.    
 * @param v    A parms_Viewer object.   
 * 
 * @return 0 on success.
 */
int parms_TimerView(parms_Timer self, parms_Viewer v)
{
  FILE *fp;

  parms_ViewerGetFP(v, &fp);
  fprintf(fp, "init_time = %f elapsed =%f\n", self->initial_time, 
	  self->elapsed_time); 

  return 0;
}

