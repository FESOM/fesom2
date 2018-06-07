SUBROUTINE nemogcmcoup_step( istp, icdate, ictime )

   IMPLICIT NONE

   ! Arguments

   ! Time step
   INTEGER, INTENT(IN) :: istp

   ! Data and time from NEMO
   INTEGER, INTENT(OUT) :: icdate, ictime

   ! Local variables
   
   ! Advance the FESOM model 1 time step

   WRITE(0,*)'Insert FESOM step here.'

   ! Compute date and time at the end of the time step.

#ifdef FESOM_TODO
   iye = ndastp / 10000
   imo = ndastp / 100 - iye * 100
   ida = MOD( ndastp, 100 )
   CALL greg2jul( 0, 0, 0, ida, imo, iye, zjul )
   zjul = zjul + ( nsec_day + 0.5_wp * rdttra(1) ) / 86400.0_wp
   CALL jul2greg( iss, imm, ihh, ida, imo, iye, zjul )
   icdate = iye * 10000 + imo * 100 + ida
   ictime = ihh * 10000 + imm * 100 + iss
#endif

END SUBROUTINE nemogcmcoup_step
   
