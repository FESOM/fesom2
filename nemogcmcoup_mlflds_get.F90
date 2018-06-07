SUBROUTINE nemogcmcoup_mlflds_get( mype, npes, icomm, &
   &                               nlev, nopoints, pgt3d, pgs3d, pgu3d, pgv3d )

   ! Interpolate sst, ice: surf T; albedo; concentration; thickness,
   ! snow thickness and currents from the ORCA grid to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   USE par_kind
   IMPLICIT NONE
   
   ! Arguments
   REAL(wp), DIMENSION(nopoints,nlev) :: pgt3d, pgs3d, pgu3d, pgv3d
   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm
   ! Number Gaussian grid points
   INTEGER, INTENT(IN) :: nopoints,nlev

   ! Local variables

   WRITE(0,*)'nemogcmcoup_mlflds_get should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_mlflds_get
   
