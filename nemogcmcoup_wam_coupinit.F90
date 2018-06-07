SUBROUTINE nemogcmcoup_wam_coupinit( mype, npes, icomm, &
   &                                 nlocpoints, nglopoints, &
   &                                 nlocmsk, ngloind, iunit )

   ! Initialize single executable coupling between WAM and NEMO
   ! This is called from WAM.

   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype,npes,icomm
   ! WAM grid information   
   ! Number of local and global points
   INTEGER, INTENT(IN) :: nlocpoints, nglopoints
   ! Integer mask and global indices
   INTEGER, DIMENSION(nlocpoints), INTENT(IN) :: nlocmsk, ngloind
   ! Unit for output in parinter_init
   INTEGER :: iunit
   
   WRITE(0,*)'Wam couplind not implemented for FESOM'
   CALL abort

END SUBROUTINE nemogcmcoup_wam_coupinit
