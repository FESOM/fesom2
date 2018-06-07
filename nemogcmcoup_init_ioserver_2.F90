SUBROUTINE nemogcmcoup_init_ioserver_2( icomm )

   ! Initialize the NEMO mppio server

   IMPLICIT NONE
   INTEGER :: icomm

   WRITE(*,*)'No mpp_ioserver'
   CALL abort
   
END SUBROUTINE nemogcmcoup_init_ioserver_2
