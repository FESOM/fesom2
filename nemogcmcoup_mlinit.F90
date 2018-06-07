SUBROUTINE nemogcmcoup_mlinit( mype, npes, icomm, &
   &                           nlev, nopoints, pdep, pmask )

   ! Get information about the vertical discretization of the ocean model
   
   ! nlevs are maximum levels on input and actual number levels on output

   USE par_kind

   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! Grid information
   INTEGER, INTENT(INOUT) :: nlev, nopoints
   REAL(wp), INTENT(OUT), DIMENSION(nlev) :: pdep
   REAL(wp), INTENT(OUT), DIMENSION(nopoints,nlev) :: pmask

   ! Local variables

   WRITE(0,*)'nemogcmcoup_mlinit should not be called when coupling to fesom.'
   CALL abort
   
END SUBROUTINE nemogcmcoup_mlinit
