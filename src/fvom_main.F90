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
IMPLICIT NONE

integer :: n, nsteps,offset,row,i
	call par_init 
        !=====================
	! Read configuration data,  
	! load the mesh and fill in 
	! auxiliary mesh arrays
	!=====================
	call setup_model          ! Read Namelists, always before clock_init
        call clock_init           ! read the clock file 
	call get_run_steps(nsteps)
	call mesh_setup
write(*,*) 'mesh_setup... complete'
   	!=====================
	! Allocate field variables 
	! and additional arrays needed for 
	! fancy advection etc.  
	!=====================
        call check_mesh_consistency
write(*,*) 'check_mesh_consistency... complete'
	call ocean_setup
write(*,*) 'ocean_setup... complete'
	call forcing_setup
	if (use_ice) then 
	  call ice_setup
          ice_steps_since_upd = ice_ave_steps-1
          ice_update=.true.
        endif
	call clock_newyear                    	! check if it is a new year
        call init_output_mean(.not. r_restart)  ! create new output files
        call init_output_restart(.not. r_restart)
	!=====================
	! Time stepping
	!=====================
write(*,*) 'start interation before the barrier...'
        call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
write(*,*) 'start interation after the barrier...'
	do n=1, nsteps
	   if (mype==0) write(*,*) 'DAY  ', n*dt/24./3600.
	   call init_output_mean(yearnew/=yearold)
	   call init_output_restart(yearnew/=yearold)
           call clock
	   call forcing_index
           call compute_vel_nodes 
	   if(use_ice) then
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
	   else
	     if(.not.toy_ocean) call update_atm_forcing_OnlyOcean(n)  
	   end if  

	   call oce_timestep(n)

	   call output (0,n)        ! save (NetCDF)
	   call restart(0,n)        ! save (NetCDF)
	end do

	if (mype==0) write(*,*) 'Run is finished, updating clock'
	call clock_finish  
	call par_ex
end program main
