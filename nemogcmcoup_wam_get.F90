SUBROUTINE nemogcmcoup_wam_get( mype, npes, icomm, &
   &                            nopoints, pwsst, pwicecov, pwicethk, &
   &                            pwucur, pwvcur, licethk )

   ! Interpolate from the ORCA grid
   ! to the WAM grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind
   IMPLICIT NONE
   
   ! Arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number WAM grid points
   INTEGER, INTENT(IN) :: nopoints
   ! Local arrays of sst, ice cover, ice thickness and currents
   REAL(wp), DIMENSION(nopoints) :: pwsst, pwicecov, pwicethk, pwucur, pwvcur
   LOGICAL :: licethk

   ! Local variables

   WRITE(0,*)'nemogcmcoup_wam_get should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_get
   
