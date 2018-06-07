SUBROUTINE nemogcmcoup_wam_update_stress( mype, npes, icomm, npoints, &
   &                                      pwutau, pwvtau, pwuv10n, pwphif,&
   &                                      cdtpro, ldebug )

   ! Update stresses in nemogcmcoup_data by parallel
   ! interpolation of the input WAM grid data
   
   USE par_kind

   IMPLICIT NONE

   ! Arguments

   ! MPI communications
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Data on the WAM grid.
   INTEGER, INTENT(IN) :: npoints
   REAL(wp), DIMENSION(npoints), INTENT(IN) :: &
      & pwutau, pwvtau, pwuv10n, pwphif
   ! Current time step
   CHARACTER(len=14), INTENT(IN) :: cdtpro
   ! Write debugging fields in netCDF
   LOGICAL, INTENT(IN) :: ldebug

   ! Local variables
   
   WRITE(0,*)'nemogcmcoup_wam_update_stress should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_update_stress

   
