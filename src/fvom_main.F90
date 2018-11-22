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

integer :: n, nsteps, offset, row, i
real(kind=WP) :: t0, t1
real(kind=real32) :: mean_rtime(9), max_rtime(9), min_rtime(9), runtime_alltimesteps
  

#ifndef __oifs
    !ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular dependencies between modules (cpl_driver and g_PARSUP)
    !OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS
    call MPI_INIT(i) 
#endif


#if defined (__oasis)
    call cpl_oasis3mct_init(MPI_COMM_FESOM)
#endif

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
    

    call clock_newyear                        ! check if it is a new year
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
    !___MIDDEPTH LEVELS NEEDS TO BE CALCULATET AT RESTART_______________________
    if (r_restart .and. use_ALE) then
        call restart_thickness_ale
    end if
    
    
    !=====================
    ! Time stepping
    !=====================
    if (mype==0) write(*,*) 'FESOM start interation before the barrier...'
    call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
    
    if (mype==0) then
       write(*,*) 'FESOM start interation after the barrier...'
       t0 = MPI_Wtime()
    endif
   
    !___MODEL TIME STEPPING LOOP________________________________________________
    do n=1, nsteps        
        mstep = n
        if (mod(n,logfile_outfreq)==0 .and. mype==0) then
            write(*,*) 'FESOM ======================================================='
!             write(*,*) 'FESOM step:',n,' day:', n*dt/24./3600.,
            write(*,*) 'FESOM step:',n,' day:', daynew,' year:',yearnew 
            write(*,*)
        end if
#if defined (__oifs) || defined (__oasis)
            seconds_til_now=INT(dt)*(n-1)
#endif
        call clock
        call compute_vel_nodes 
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
    
    !___FINISH MODEL RUN________________________________________________________

    call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
    if (mype==0) then
       t1 = MPI_Wtime()
       runtime_alltimesteps = real(t1-t0,real32)
       write(*,*) 'FESOM Run is finished, updating clock'
    endif
    
    mean_rtime(1) = rtime_oce         
    mean_rtime(2) = rtime_oce_mixpres 
    mean_rtime(3) = rtime_oce_dyn     
    mean_rtime(4) = rtime_oce_dynssh  
    mean_rtime(5) = rtime_oce_solvessh
    mean_rtime(6) = rtime_oce_GMRedi  
    mean_rtime(7) = rtime_oce_solvetra
    mean_rtime(8) = rtime_ice         
    mean_rtime(9) = rtime_tot  

    max_rtime(1:9) = mean_rtime(1:9)
    min_rtime(1:9) = mean_rtime(1:9)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 9, MPI_REAL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    mean_rtime(1:9) = mean_rtime(1:9) / real(npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  9, MPI_REAL, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  9, MPI_REAL, MPI_MIN, MPI_COMM_FESOM, MPIerr)

    if (mype==0) then
        write(*,*) '___MODEL RUNTIME mean, min, max per task [seconds]________________________'
        write(*,*) '  runtime ocean:',mean_rtime(1), min_rtime(1), max_rtime(1)
        write(*,*) '    > runtime oce. mix,pres. :',mean_rtime(2), min_rtime(2), max_rtime(2)
        write(*,*) '    > runtime oce. dyn. u,v,w:',mean_rtime(3), min_rtime(3), max_rtime(3)
        write(*,*) '    > runtime oce. dyn. ssh  :',mean_rtime(4), min_rtime(4), max_rtime(4)
        write(*,*) '        > runtime oce. solve ssh:',mean_rtime(5), min_rtime(5), max_rtime(5)
        write(*,*) '    > runtime oce. GM/Redi   :',mean_rtime(6), min_rtime(6), max_rtime(6)
        write(*,*) '    > runtime oce. tracer    :',mean_rtime(7), min_rtime(7), max_rtime(7)
        write(*,*) '  runtime ice  :',mean_rtime(8), min_rtime(8), max_rtime(8)
        write(*,*) '  runtime total:',mean_rtime(9), min_rtime(9), max_rtime(9)
        write(*,*)
        write(*,*) '============================================'
        write(*,*) '=========== BENCHMARK RUNTIME =============='
        write(*,*) '    Number of cores : ',npes
        write(*,*) '    Runtime for all timesteps : ',runtime_alltimesteps,' sec'
        write(*,*) '============================================'
        write(*,*)
    end if    
!   call clock_finish  
    call par_ex
end program main
    
