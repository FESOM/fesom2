!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    

program main

  use g_PARSUP, only: mype
  integer :: nsteps

  call main_initialize(nsteps)
  if (mype==0) write(*,*) 'Initialization complete.'

  call main_timestepping(nsteps)
  if (mype==0) write(*,*) 'Timestepping complete...'

  call main_finalize
  if (mype==0) write(*,*) 'Finalization complete...'

end program main



!=============================================================================!
!                      Initialization
!=============================================================================!

subroutine main_initialize(nsteps)
  ! Split main into three major parts
  ! Coded by Thomas Rackow, 2018
  !----------------------------------
  USE g_PARSUP, only: mype, par_init
  USE i_PARAM,  only: ice_ave_steps, whichEVP
  use i_ARRAYS, only: ice_steps_since_upd, ice_update
  use g_clock,  only: clock_init, clock_newyear
  use g_config, only: use_ice, r_restart, use_ALE
  use io_RESTART,   only: restart
  use io_mesh_info, only: write_mesh_info
  use diagnostics,  only: compute_diagnostics
#if defined (__oasis)
  use cpl_driver
#endif
  IMPLICIT NONE

  integer :: ierr
  integer, INTENT(OUT) :: nsteps

#ifndef __oifs
  	!ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular 		!dependencies between modules (cpl_driver and g_PARSUP)
	!OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS
        call MPI_INIT(ierr) 
#endif


#if defined (__oasis)
        call cpl_oasis3mct_init(MPI_COMM_FESOM)
#endif

	! sets npes and mype
	call par_init 
	if (mype==0) write(*,*) '!=============================================================================!'
	if (mype==0) write(*,*) '!                             Welcome to the					'
	if (mype==0) write(*,*) '!                 Finite Volume Sea-ice Ocean Model (FESOM2)                   '
	if (mype==0) write(*,*) '!=============================================================================!'

	!=====================
	! Read configuration data,  
	! load the mesh and fill in 
	! auxiliary mesh arrays
	!=====================
	call setup_model          ! Read Namelists, always before clock_init
	call clock_init           ! read the clock file 
	call get_run_steps(nsteps)
	call mesh_setup
	if (mype==0) write(*,*) 'FESOM mesh_setup... complete'
	
	!=====================
	! Allocate field variables 
	! and additional arrays needed for 
	! fancy advection etc.  
	!=====================
	call check_mesh_consistency
	call ocean_setup
	if (mype==0) write(*,*) 'FESOM ocean_setup... complete'
	call forcing_setup
	if (use_ice) then 
		call ice_setup
		ice_steps_since_upd = ice_ave_steps-1
		ice_update=.true.
                if (mype==0) write(*,*) 'EVP scheme option=', whichEVP
	endif
        call compute_diagnostics(0) ! allocate arrays for diagnostic
#if defined (__oasis)
	call cpl_oasis3mct_define_unstr
	if(mype==0)  write(*,*) 'FESOM ---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv
#endif
	

	call clock_newyear                    	! check if it is a new year
	!___CREATE NEW RESTART FILE IF APPLICABLE___________________________________
	! The interface to the restart module is made via call restart !
	! The inputs are: istep, l_write, l_create
	! if istep is not zero it will be decided whether restart shall be written
	! if l_write  is TRUE the restart will be forced
	! if l_read the restart will be read
	! as an example, for reading restart one does: call restart(0, .false., .false., .true.)
	call restart(0, .false., r_restart) ! istep, l_write, l_read
        ! store grid information into netcdf file
	if (.not. r_restart) call write_mesh_info

	!___IF RESTART WITH ZLEVEL OR ZSTAR IS DONE, ALSO THE ACTUAL LEVELS AND ____
	!___MIDDEPTH LEVELS NEEDS TO BE CALCULATED AT RESTART_______________________
	if (r_restart .and. use_ALE) then
		call restart_thickness_ale
	end if

end subroutine main_initialize



!=============================================================================!
!                      Timestepping
!=============================================================================!

subroutine main_timestepping(nsteps)
  ! Split main into three major parts
  ! Coded by Thomas Rackow, 2018
  !----------------------------------
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
  use io_RESTART
  use io_MEANDATA
  use io_mesh_info
  use diagnostics
#if defined (__oasis)
  use cpl_driver
#endif
  IMPLICIT NONE

  integer :: n
  integer, INTENT(IN) :: nsteps 

	!=====================
	! Time stepping
	!=====================
	if (mype==0) write(*,*) 'FESOM start integration before the barrier...'
		call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
	if (mype==0) write(*,*) 'FESOM start integration after the barrier...'
	
	
	!___MODEL TIME STEPPING LOOP________________________________________________
	do n=1, nsteps
		
		mstep = n
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,*) 'FESOM ======================================================='
! 			write(*,*) 'FESOM step:',n,' day:', n*dt/24./3600.,
			write(*,*) 'FESOM step:',n,' day:', daynew,' year:',yearnew 
			write(*,*)
		end if
#if defined (__oifs) || defined (__oasis)
        	seconds_til_now=INT(dt)*(n-1)
#endif
		call clock
		call forcing_index
		call compute_vel_nodes 
		
! 		eta_n=alpha*hbar+(1.0_WP-alpha)*hbar_old !PS
		
		!___model sea-ice step__________________________________________________
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
			
			call oce_fluxes_mom ! momentum only
                        call oce_fluxes
		end if  
		
		

		!___model ocean step____________________________________________________
		call oce_timestep_ale(n)
                call compute_diagnostics(1)
		!___prepare output______________________________________________________
		call output (n)
		call restart(n, .false., .false.)
	end do
end subroutine main_timestepping



!=============================================================================!
!                      Finalization
!=============================================================================!

subroutine main_finalize
  ! Split main into three major parts
  ! Coded by Thomas Rackow, 2018
  !----------------------------------

  USE g_PARSUP, only: mype, par_ex
  use g_clock, only: clock_finish

  IMPLICIT NONE

  !___FINISH MODEL RUN________________________________________________________
  if (mype==0) write(*,*) 'FESOM run is finished, updating clock'
  call clock_finish  
  call par_ex

end subroutine main_finalize
