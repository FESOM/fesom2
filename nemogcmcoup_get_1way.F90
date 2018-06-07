SUBROUTINE nemogcmcoup_get_1way( mype, npes, icomm )

   ! Interpolate sst, ice and currents from the ORCA grid
   ! to the Gaussian grid. 
   
   ! This routine can be called at any point in time since it does
   ! the necessary message passing in parinter_fld. 

   IMPLICIT NONE

   
   ! Arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: mype, npes, icomm

   ! Local variables

   WRITE(0,*)'nemogcmcoup_get_1way should not be called when coupling to fesom.'
   CALL abort

END SUBROUTINE nemogcmcoup_get_1way
   
