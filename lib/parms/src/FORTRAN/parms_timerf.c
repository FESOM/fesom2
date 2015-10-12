#include "parms_timer.h"

#if defined(FORTRAN_CAPS)
#define parms_timercreate_          PARMS_TIMERCREATE      
#define parms_timerfree_ 	    PARMS_TIMERFREE        
#define parms_timerget_ 	    PARMS_TIMERGET         
#define parms_timerpause_ 	    PARMS_TIMERPAUSE       
#define parms_timerreset_ 	    PARMS_TIMERRESET       
#define parms_timerresetdelay_ 	    PARMS_TIMERRESETDELAY  
#define parms_timerrestart_ 	    PARMS_TIMERRESTART     
#define parms_timerview_            PARMS_TIMERVIEW        
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define parms_timercreate_          parms_timercreate__      
#define parms_timerfree_ 	    parms_timerfree__        
#define parms_timerget_ 	    parms_timerget__         
#define parms_timerpause_ 	    parms_timerpause__       
#define parms_timerreset_ 	    parms_timerreset__       
#define parms_timerresetdelay_ 	    parms_timerresetdelay__  
#define parms_timerrestart_ 	    parms_timerrestart__     
#define parms_timerview_            parms_timerview__        
#elif !defined(FORTRAN_UNDERSCORE)
#define parms_timercreate_          parms_timercreate      
#define parms_timerfree_ 	    parms_timerfree        
#define parms_timerget_ 	    parms_timerget         
#define parms_timerpause_ 	    parms_timerpause       
#define parms_timerreset_ 	    parms_timerreset       
#define parms_timerresetdelay_ 	    parms_timerresetdelay  
#define parms_timerrestart_ 	    parms_timerrestart     
#define parms_timerview_            parms_timerview        
#endif

void parms_timercreate_(parms_Timer *self, int *ierr)
{
  parms_TimerCreate(self);
  *ierr = 0;
}

void parms_timerfree_(parms_Timer *self, int *ierr)
{
  *ierr = parms_TimerFree(self);
}

void parms_timerget_(parms_Timer *self, double *t, int *ierr)
{
  *t = parms_TimerGet(*self);
  *ierr = 0;
}

void parms_timerpause_(parms_Timer *self, int *ierr)
{
  *ierr = parms_TimerPause(*self);
}

void parms_timerreset_(parms_Timer *self, int *ierr)    
{
  *ierr = parms_TimerReset(*self);
}

void parms_timerresetdelay_(parms_Timer *self, double *delay, int
			    *ierr)
{
  parms_TimerResetDelay(*self, *delay);
  *ierr = 0;
}

void parms_timerrestart_(parms_Timer *self, int *ierr)
{
  *ierr = parms_TimerRestart(*self);
}

void parms_timerview_(parms_Timer *self, parms_Viewer *v, int *ierr)
{
  *ierr = parms_TimerView(*self, *v);
}
