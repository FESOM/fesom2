SUBROUTINE nemogcmcoup_update_add( mype, npes, icomm, &
   &                               npoints, pgsst, pgtsk, kt, ldebug )

   ! Update addetiona in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Input on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wp), DIMENSION(npoints), intent(IN) :: &
      & pgsst, pgtsk
   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables

   WRITE(0,*)'nemogcmcoup_update_add should not be called when coupling to fesom.'
   CALL abort
   

END SUBROUTINE nemogcmcoup_update_add

   
