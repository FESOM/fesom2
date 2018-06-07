SUBROUTINE nemogcmcoup_get( mype, npes, icomm, &
   &                        nopoints, pgsst, pgice, pgucur, pgvcur )

   ! Interpolate sst, ice and currents from the ORCA grid
   ! to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind

   IMPLICIT NONE

   
   ! Arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints
   ! Local arrays of sst, ice and currents
   REAL(wp), DIMENSION(nopoints) :: pgsst, pgice, pgucur, pgvcur

   ! Local variables

   WRITE(0,*)'nemogcmcoup_get should not be called with FESOM'
   CALL abort

END SUBROUTINE nemogcmcoup_get
   
