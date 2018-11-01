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
use io_RESTART
use io_MEANDATA
use io_mesh_info
use diagnostics
#if defined (__oasis)
use cpl_driver
#endif
IMPLICIT NONE

integer :: n, nsteps, offset, row, i
real(kind=WP) :: mrtime_ice=0.0,mrtime_oce=0.0,mrtime_tot=0.0
real(kind=WP) :: mrtime_oce_dyn=0.0, mrtime_oce_dynssh=0.0, mrtime_oce_solvessh=0.0 
real(kind=WP) :: mrtime_oce_solvetra=0.0, mrtime_oce_GMRedi=0.0, mrtime_oce_mixpres=0.0
  

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
    if (mype==0) write(*,*) 'FESOM start interation after the barrier...'
   
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
        call forcing_index
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
    if (mype==0) write(*,*) 'FESOM Run is finished, updating clock'
    
    ! average ocean, ice and total runtime over all cpus
    call MPI_AllREDUCE(rtime_oce         , mrtime_oce         , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_oce_mixpres , mrtime_oce_mixpres , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_oce_dyn     , mrtime_oce_dyn     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_oce_dynssh  , mrtime_oce_dynssh  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_oce_solvessh, mrtime_oce_solvessh, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_oce_GMRedi  , mrtime_oce_GMRedi  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_oce_solvetra, mrtime_oce_solvetra, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_ice         , mrtime_ice         , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(rtime_tot         , mrtime_tot         , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    if (mype==0) then
        write(*,*) '___MODEL RUNTIME [seconds]_____________________________'
        write(*,*) '    runtime ocean : ',mrtime_oce/npes, ' sec'
        write(*,*) '      > runtime oce. mix,pres...: ',mrtime_oce_mixpres/npes, ' sec'
        write(*,*) '      > runtime oce. dyn. u,v,w : ',mrtime_oce_dyn/npes, ' sec'
        write(*,*) '      > runtime oce. dyn. ssh   : ',mrtime_oce_dynssh/npes, ' sec'
        write(*,*) '          > runtime oce. solve ssh  : ',mrtime_oce_solvessh/npes, ' sec'
        write(*,*) '      > runtime oce. GM/Redi    : ',mrtime_oce_GMRedi/npes, ' sec'
        write(*,*) '      > runtime oce. solve tacer: ',mrtime_oce_solvetra/npes, ' sec'
        write(*,*) '    runtime ice   : ',mrtime_ice/npes, ' sec'
        write(*,*) '    runtime total : ',mrtime_tot/npes, ' sec'
    end if     
!   call clock_finish  
    call par_ex
end program main
    
