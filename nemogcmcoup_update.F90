SUBROUTINE nemogcmcoup_update( mype, npes, icomm, &
   &                           npoints, pgutau, pgvtau, &
   &                           pgqsr, pgqns, pgemp, kt, ldebug )

   ! Update fluxes in nemogcmcoup_data by parallel
   ! interpolation of the input gaussian grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Fluxes on the Gaussian grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wp), DIMENSION(npoints), intent(IN) :: &
      & pgutau, pgvtau, pgqsr, pgqns, pgemp
   ! Current time step
   INTEGER, INTENT(in) :: kt
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables

   WRITE(0,*)'nemogcmcoup_update should be called with with.'
   CALL abort
   
END SUBROUTINE nemogcmcoup_update

   
