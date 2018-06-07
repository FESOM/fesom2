SUBROUTINE nemogcmcoup_exflds_get( mype, npes, icomm, &
   &                               nopoints, pgssh, pgmld, pg20d, pgsss, &
   &                               pgtem300, pgsal300 )

   ! Interpolate sst, ice: surf T; albedo; concentration; thickness,
   ! snow thickness and currents from the ORCA grid to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind
   IMPLICIT NONE
   
   ! Arguments
   REAL(wp), DIMENSION(nopoints) :: pgssh, pgmld, pg20d, pgsss, &
      & pgtem300, pgsal300
   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints

   ! Local variables

   WRITE(0,*)'nemogcmcoup_exflds_get should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_exflds_get
   
