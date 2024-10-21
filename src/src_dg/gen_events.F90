!
!--------------------------------------------------------------------------------------------
!
subroutine annual_event(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if ((daynew == ndpyr) .and. (timenew==86400.)) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine annual_event
!
!--------------------------------------------------------------------------------------------
!
subroutine monthly_event(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (day_in_month==num_day_in_month(fleapyear,month) .and. &
       timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  end if

end subroutine monthly_event
!
!--------------------------------------------------------------------------------------------
!
subroutine daily_event(do_output, N)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical             :: do_output
  integer, intent(in) :: N
  if (mod(daynew, N)==0 .and. timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine daily_event
!
!--------------------------------------------------------------------------------------------
!
subroutine hourly_event(do_output, N)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical             :: do_output
  integer, intent(in) :: N

  if (mod(timenew, 3600.*N)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine hourly_event
!
!--------------------------------------------------------------------------------------------
!
subroutine step_event(do_output,istep, N)
  !decides whether it's time to do output
  use g_config
  implicit none

  logical             :: do_output
  integer             :: istep
  integer, intent(in) :: N

  if (mod(istep, N)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine step_event
!
!--------------------------------------------------------------------------------------------
!
subroutine handle_err(errcode)
  use g_parsup
  implicit none
  
#include "netcdf.inc" 
  
  integer errcode
  
  write(*,*) 'Error: ', nf_strerror(errcode)
  call par_ex(1)
  stop
end subroutine handle_err
!
!--------------------------------------------------------------------------------------------
!
