!=============================================================================!
!
!                 Finite Volume Sea-ice Ocean Model
!
!=============================================================================!
!                      The main driving routine
!=============================================================================!    

program main
USE MOD_MESH
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
use mo_tidal
use fesom_version_info_module


! Define icepack module
#if defined (__icepack)
use icedrv_main,          only: set_icepack, init_icepack, alloc_icepack
#endif

#if defined (__oasis)
use cpl_driver
#endif

#if defined(__recom)
  use REcoM_GloVar
  use recom_config
  use recom_diag
! for diagnosing si mass balance:
  use g_support                                                                                               
#endif

IMPLICIT NONE

integer :: n, nsteps, offset, row, i, provided
real(kind=WP)     :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t0_ice, t1_ice, t0_frc, t1_frc
real(kind=WP)     :: rtime_fullice,    rtime_write_restart, rtime_write_means, rtime_compute_diag, rtime_read_forcing
real(kind=real32) :: rtime_setup_mesh, rtime_setup_ocean, rtime_setup_forcing 
real(kind=real32) :: rtime_setup_ice,  rtime_setup_other, rtime_setup_restart
real(kind=real32) :: mean_rtime(15), max_rtime(15), min_rtime(15)
real(kind=real32) :: runtime_alltimesteps
#if defined(__recom)
real(kind=WP),  save,  target                 :: intDSi
real(kind=WP),  save,  target                 :: intDiaSi
real(kind=WP),  save,  target                 :: intDetSi
real(kind=WP),  save,  target                 :: intBenSi
real(kind=WP),  save,  target                 :: sumSi1, sumSi2
#endif

type(t_mesh),             target, save :: mesh

#ifndef __oifs
    !ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular dependencies between modules (cpl_driver and g_PARSUP)
    !OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS
    call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, i)
#endif
    

#if defined (__oasis)
    call cpl_oasis3mct_init(MPI_COMM_FESOM)
#endif
    t1 = MPI_Wtime()

    call par_init 
    if(mype==0) then
        write(*,*)
        print *,"FESOM2 git SHA: "//fesom_git_sha()
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> FESOM BUILDS UP MODEL CONFIGURATION                    '//achar(27)//'[0m'
    end if
    !=====================
    ! Read configuration data,  
    ! load the mesh and fill in 
    ! auxiliary mesh arrays
    !=====================
    call setup_model          ! Read Namelists, always before clock_init
    call clock_init           ! read the clock file 
    call get_run_steps(nsteps)
    call mesh_setup(mesh)

    if (mype==0) write(*,*) 'FESOM mesh_setup... complete'
    
    !=====================
    ! Allocate field variables 
    ! and additional arrays needed for 
    ! fancy advection etc.  
    !=====================
    call check_mesh_consistency(mesh)
    if (mype==0) t2=MPI_Wtime()
    call ocean_setup(mesh)
    if (mype==0) then
       write(*,*) 'FESOM ocean_setup... complete'
       t3=MPI_Wtime()
    endif
#if defined (__recom)
       call recom_init(mesh)
       if (mype==0) write(*,*) 'RECOM recom_init... complete'
#endif    
    call forcing_setup(mesh)
    if (mype==0) t4=MPI_Wtime()
    if (use_ice) then 
        call ice_setup(mesh)
        ice_steps_since_upd = ice_ave_steps-1
        ice_update=.true.
        if (mype==0) write(*,*) 'EVP scheme option=', whichEVP
    endif
#if defined(__recom)
        call compute_recom_diagnostics(0, mesh) ! allocate arrays for recom diagnostic
#endif

    if (mype==0) t5=MPI_Wtime()
    call compute_diagnostics(0, mesh) ! allocate arrays for diagnostic
#if defined (__oasis)
    call cpl_oasis3mct_define_unstr(mesh)
    if(mype==0)  write(*,*) 'FESOM ---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv
#endif

#if defined (__icepack)
    !=====================
    ! Setup icepack
    !=====================
    if (mype==0) write(*,*) 'Icepack: reading namelists from namelist.icepack'
    call set_icepack
    call alloc_icepack
    call init_icepack(mesh)
    if (mype==0) write(*,*) 'Icepack: setup complete'
#endif
    
    call clock_newyear                        ! check if it is a new year
    if (mype==0) t6=MPI_Wtime()
    !___CREATE NEW RESTART FILE IF APPLICABLE___________________________________
    ! The interface to the restart module is made via call restart !
    ! The inputs are: istep, l_write, l_create
    ! if istep is not zero it will be decided whether restart shall be written
    ! if l_write  is TRUE the restart will be forced
    ! if l_read the restart will be read
    ! as an example, for reading restart one does: call restart(0, .false., .false., .true.)
    call restart(0, .false., r_restart, mesh) ! istep, l_write, l_read
    if (mype==0) t7=MPI_Wtime()
    
    ! store grid information into netcdf file
    if (.not. r_restart) call write_mesh_info(mesh)

    !___IF RESTART WITH ZLEVEL OR ZSTAR IS DONE, ALSO THE ACTUAL LEVELS AND ____
    !___MIDDEPTH LEVELS NEEDS TO BE CALCULATET AT RESTART_______________________
    if (r_restart) then
        call restart_thickness_ale(mesh)
    end if

    if (mype==0) then
       t8=MPI_Wtime()
    
       rtime_setup_mesh    = real( t2 - t1              ,real32)
       rtime_setup_ocean   = real( t3 - t2              ,real32)
       rtime_setup_forcing = real( t4 - t3              ,real32)
       rtime_setup_ice     = real( t5 - t4              ,real32)
       rtime_setup_restart = real( t7 - t6              ,real32)
       rtime_setup_other   = real((t8 - t7) + (t6 - t5) ,real32)

       write(*,*) '=========================================='
       write(*,*) 'MODEL SETUP took on mype=0 [seconds]      '
       write(*,*) 'runtime setup total      ',real(t8-t1,real32)      
       write(*,*) ' > runtime setup mesh    ',rtime_setup_mesh   
       write(*,*) ' > runtime setup ocean   ',rtime_setup_ocean  
       write(*,*) ' > runtime setup forcing ',rtime_setup_forcing
       write(*,*) ' > runtime setup ice     ',rtime_setup_ice    
       write(*,*) ' > runtime setup restart ',rtime_setup_restart
       write(*,*) ' > runtime setup other   ',rtime_setup_other 
        write(*,*) '============================================' 
    endif

    !=====================
    ! Time stepping
    !=====================

! Initialize timers
    rtime_fullice       = 0._WP
    rtime_write_restart = 0._WP
    rtime_write_means   = 0._WP
    rtime_compute_diag  = 0._WP
    rtime_read_forcing  = 0._WP

    if (mype==0) write(*,*) 'FESOM start iteration before the barrier...'
    call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
    
    if (mype==0) then
       write(*,*) 'FESOM start iteration after the barrier...'
       t0 = MPI_Wtime()
    endif
    if(mype==0) then
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> FESOM STARTS TIME LOOP                                 '//achar(27)//'[0m'
    end if
    !___MODEL TIME STEPPING LOOP________________________________________________
    if (use_global_tides) then
       call foreph_ini(yearnew, month)
    end if

    do n=1, nsteps        
        if (use_global_tides) then
           call foreph(mesh)
        end if
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
        
        !___compute horizontal velocity on nodes (originaly on elements)________
        call compute_vel_nodes(mesh)
        
        !___model sea-ice step__________________________________________________
        t1 = MPI_Wtime()
        if(use_ice) then
            !___compute fluxes from ocean to ice________________________________
            if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call ocean2ice(n)'//achar(27)//'[0m'
            call ocean2ice(mesh)
            
            !___compute update of atmospheric forcing____________________________
            if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing(n)'//achar(27)//'[0m'
            t0_frc = MPI_Wtime()
            call update_atm_forcing(n, mesh)
            t1_frc = MPI_Wtime()            
            !___compute ice step________________________________________________
            if (ice_steps_since_upd>=ice_ave_steps-1) then
                ice_update=.true.
                ice_steps_since_upd = 0
            else
                ice_update=.false.
                ice_steps_since_upd=ice_steps_since_upd+1
            endif
            if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
            if (ice_update) call ice_timestep(n, mesh)  
            !___compute fluxes to the ocean: heat, freshwater, momentum_________
            if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call oce_fluxes_mom...'//achar(27)//'[0m'
            call oce_fluxes_mom(mesh) ! momentum only
            call oce_fluxes(mesh)
        end if
        call before_oce_step(mesh) ! prepare the things if required
#if defined (__recom)
        if (use_REcoM) then
           call recom(mesh)



       ! for silicate mass balance:
!        if (mype==0) print *, '2) si tracer global integral in fvom_main before oce_timestep_ale:'
!        call integrate_nod(tr_arr(:,:,isi+2), intDSi, mesh)
!        call integrate_nod(tr_arr(:,:,idiasi+2), intDiaSi, mesh)
!        call integrate_nod(tr_arr(:,:,idetsi+2), intDetSi, mesh)
!        call integrate_nod(Benthos(:,3), intBenSi, mesh)
        !if (mype==0) print *, 'intDSi: ', intDSi
        !if (mype==0) print *, 'intDiaSi: ', intDiaSi
        !if (mype==0) print *, 'intDetSi: ', intDetSi
        !if (mype==0) print *, 'intBenSi: ', intBenSi
!        sumSi1 = intDSi + intDiaSi + intDetSi + intBenSi
!        if (mype==0 .and. (sumSi1-sumSi8)<0) print *, 'sumSi8, sumSi1, fvom before oce_step: ', sumSi8, sumSi1
!        if (mype==0) print *, 'sumSi1: ', sumSi1
!        if (mype==0) print *, ''
!        if (mype==0) print *, ''

           if (mype==0 .and. n==1)  print *, achar(27)//'[46;1m'//'     --> call RECOM         '//achar(27)//'[0m'
           call compute_recom_diagnostics(1, mesh)
        end if
#endif

        t2 = MPI_Wtime()
        
        !___model ocean step____________________________________________________
        if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call oce_timestep_ale'//achar(27)//'[0m'
        call oce_timestep_ale(n, mesh)
        t3 = MPI_Wtime()
        !___compute energy diagnostics..._______________________________________
        if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics(1)'//achar(27)//'[0m'
        call compute_diagnostics(1, mesh)
        t4 = MPI_Wtime()
        !___prepare output______________________________________________________
        if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call output (n)'//achar(27)//'[0m'
        call output (n, mesh)
        t5 = MPI_Wtime()
        call restart(n, .false., .false., mesh)
        t6 = MPI_Wtime()
        
        rtime_fullice       = rtime_fullice       + t2 - t1
        rtime_compute_diag  = rtime_compute_diag  + t4 - t3
        rtime_write_means   = rtime_write_means   + t5 - t4   
        rtime_write_restart = rtime_write_restart + t6 - t5
        rtime_read_forcing  = rtime_read_forcing  + t1_frc - t0_frc
    end do
    
    call finalize_output()
    
    !___FINISH MODEL RUN________________________________________________________

    call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
    if (mype==0) then
       t1 = MPI_Wtime()
       runtime_alltimesteps = real(t1-t0,real32)
       write(*,*) 'FESOM Run is finished, updating clock'
    endif
    
    mean_rtime(1)  = rtime_oce         
    mean_rtime(2)  = rtime_oce_mixpres 
    mean_rtime(3)  = rtime_oce_dyn     
    mean_rtime(4)  = rtime_oce_dynssh  
    mean_rtime(5)  = rtime_oce_solvessh
    mean_rtime(6)  = rtime_oce_GMRedi  
    mean_rtime(7)  = rtime_oce_solvetra
    mean_rtime(8)  = rtime_ice         
    mean_rtime(9)  = rtime_tot  
    mean_rtime(10) = rtime_fullice - rtime_read_forcing 
    mean_rtime(11) = rtime_compute_diag
    mean_rtime(12) = rtime_write_means
    mean_rtime(13) = rtime_write_restart
    mean_rtime(14) = rtime_read_forcing   
    
    max_rtime(1:14) = mean_rtime(1:14)
    min_rtime(1:14) = mean_rtime(1:14)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 14, MPI_REAL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    mean_rtime(1:14) = mean_rtime(1:14) / real(npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  14, MPI_REAL, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  14, MPI_REAL, MPI_MIN, MPI_COMM_FESOM, MPIerr)

    if (mype==0) then
        write(*,*) '___MODEL RUNTIME mean, min, max per task [seconds]________________________'
        write(*,*) '  runtime ocean:',mean_rtime(1), min_rtime(1), max_rtime(1)
        write(*,*) '    > runtime oce. mix,pres. :',mean_rtime(2), min_rtime(2), max_rtime(2)
        write(*,*) '    > runtime oce. dyn. u,v,w:',mean_rtime(3), min_rtime(3), max_rtime(3)
        write(*,*) '    > runtime oce. dyn. ssh  :',mean_rtime(4), min_rtime(4), max_rtime(4)
        write(*,*) '        > runtime oce. solve ssh:',mean_rtime(5), min_rtime(5), max_rtime(5)
        write(*,*) '    > runtime oce. GM/Redi   :',mean_rtime(6), min_rtime(6), max_rtime(6)
        write(*,*) '    > runtime oce. tracer    :',mean_rtime(7), min_rtime(7), max_rtime(7)
        write(*,*) '  runtime ice  :',mean_rtime(10), min_rtime(10), max_rtime(10)
        write(*,*) '    > runtime ice step :',mean_rtime(8), min_rtime(8), max_rtime(8)
        write(*,*) '  runtime diag:   ', mean_rtime(11), min_rtime(11), max_rtime(11)
        write(*,*) '  runtime output: ', mean_rtime(12), min_rtime(12), max_rtime(12)
        write(*,*) '  runtime restart:', mean_rtime(13), min_rtime(13), max_rtime(13)
        write(*,*) '  runtime forcing:', mean_rtime(14), min_rtime(14), max_rtime(14)
        write(*,*) '  runtime total (ice+oce):',mean_rtime(9), min_rtime(9), max_rtime(9)
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
    
