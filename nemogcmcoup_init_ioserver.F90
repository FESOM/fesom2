SUBROUTINE nemogcmcoup_init_ioserver( icomm, lnemoioserver )

   ! Initialize the NEMO mppio server

   IMPLICIT NONE
   INTEGER :: icomm
   LOGICAL :: lnemoioserver

   WRITE(*,*)'No mpp_ioserver'
   CALL abort
   
END SUBROUTINE nemogcmcoup_init_ioserver
