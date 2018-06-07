SUBROUTINE nemogcmcoup_init( icomm, inidate, initime, itini, itend, zstp, &
   & lwaveonly, iatmunit, lwrite )

   ! Initialize the NEMO model for single executable coupling 

   USE par_kind

   IMPLICIT NONE

   ! Input arguments

   ! Message passing information
   INTEGER, INTENT(IN) :: icomm
   ! Initial date, time, initial timestep and final time step
   INTEGER, INTENT(OUT) ::  inidate, initime, itini, itend
   ! Length of the time step
   REAL(wp), INTENT(OUT) :: zstp
   ! Coupling to waves only
   LOGICAL, INTENT(IN) :: lwaveonly
   ! Logfile unit (used if >=0)
   INTEGER :: iatmunit
   ! Write to this unit
   LOGICAL :: lwrite

   WRITE(0,*)'Insert FESOM init here.'
   CALL abort

   ! Set information for the caller

#ifdef FESOM_TODO
   inidate = nn_date0
   initime = nn_time0*3600
   itini   = nit000
   itend   = nn_itend
   zstp    = rdttra(1)
#endif

END SUBROUTINE nemogcmcoup_init
