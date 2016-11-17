!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    

program main
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
USE i_PARAM
use i_ARRAYS
use g_clock
use g_config
use g_forcing_index
use g_comm_auto
use g_forcing_arrays
USE TIME_MEASURES, ONLY : timer_init, timer, timer_finalize
IMPLICIT NONE

integer :: n, nsteps,offset,row,i
	call par_init 

 ! Initialize time measurements
 CALL TIMER_INIT
 CALL TIMER( 1, 'start', NAME = 'Total' )

        !=====================
	! Read configuration data,  
	! load the mesh and fill in 
	! auxiliary mesh arrays
	!=====================
 CALL TIMER( 2, 'start', NAME = 'Model initialization' )

	call setup_model          ! Read Namelists, always before clock_init
        call clock_init           ! read the clock file 
	call get_run_steps(nsteps)
	call mesh_setup
if (mype==0) write(*,*) mype, 'mesh_setup... complete'
   	!=====================
	! Allocate field variables 
	! and additional arrays needed for 
	! fancy advection etc.  
	!=====================
        call check_mesh_consistency
if (mype==0) write(*,*) 'check_mesh_consistency... complete'
	call ocean_setup
if (mype==0) write(*,*) 'ocean_setup... complete'
	call forcing_setup
	if (use_ice) then 
	  call ice_setup
          ice_steps_since_upd = ice_ave_steps-1
          ice_update=.true.
        endif
	call clock_newyear                    	! check if it is a new year
        call init_output_mean(.not. r_restart)  ! create new output files
        call init_output_restart(.not. r_restart)

  CALL TIMER( 2, 'stop', NAME = 'Model initialization' )
        
  CALL TIMER( 3, 'start', NAME = 'Time stepping' )

	!=====================
	! Time stepping
	!=====================
if (mype==0) write(*,*) 'start integration before the barrier...'
        call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
if (mype==0) write(*,*) 'start integration after the barrier...'
	do n=1, nsteps
	   if (mype==0) write(*,*) 'DAY  ', n*dt/24./3600.
	   call init_output_mean(yearnew/=yearold)
	   call init_output_restart(yearnew/=yearold)
           call clock
	   call forcing_index
           call compute_vel_nodes 
	   if(use_ice) then
 CALL TIMER( 4, 'start', NAME = 'Ice stepping' )
	     call ocean2ice
             call update_atm_forcing(n)
             if (ice_steps_since_upd>=ice_ave_steps-1) then
              ice_update=.true.
              ice_steps_since_upd = 0
             else
              ice_update=.false.
              ice_steps_since_upd=ice_steps_since_upd+1
             endif

             if (ice_update) call ice_timestep(n)
             call ice2ocean
 CALL TIMER( 4, 'stop', NAME = 'Ice stepping' )
	   else
	     if(.not.toy_ocean) call update_atm_forcing_OnlyOcean(n)  
	   end if  

 CALL TIMER( 5, 'start', NAME = 'Ocean stepping' )
	   call oce_timestep(n)
 CALL TIMER( 5, 'stop', NAME = 'Ocean stepping' )

	   call output (0,n)        ! save (NetCDF)
	   call restart(0,n)        ! save (NetCDF)
	end do

 CALL TIMER( 3, 'stop', NAME = 'Time stepping' )

	if (mype==0) write(*,*) 'Run is finished, updating clock'
	call clock_finish  

 CALL checksums

 ! Additional timing informations
 CALL timer( 1, 'stop', NAME = 'Total' )
 CALL timer_finalize( 'FESOM2 time measurements' )

	call par_ex
end program main
