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
use iceberg_params
use mo_tidal
use fesom_version_info_module

! Define icepack module
#if defined (__icepack)
use icedrv_main,          only: set_icepack, init_icepack, alloc_icepack
#endif

#if defined (__oasis)
use cpl_driver
#endif

! kh 02.02.21
use omp_lib

IMPLICIT NONE

! kh 01.03.21
integer :: n, n_ib, nsub, nsteps, offset, row, i, provided

real(kind=WP)     :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t0_ice, t1_ice, t0_frc, t1_frc, t1_icb, t1b_icb, t2_icb, t2b_icb, t3_icb, t4_icb

! kh 09.02.21 
real(kind=WP)     :: t1_1st_section, t2_1st_section, t1_2nd_section, t2_2nd_section

! kh 24.02.21
real(kind=WP)     :: t1_par_sections, t2_par_sections

real(kind=WP)     :: rtime_fullice,    rtime_write_restart, rtime_write_means, rtime_compute_diag, rtime_icb_calc, rtime_icb_write, rtime_read_forcing

! kh 09.02.21
real(kind=WP)     :: time_1st_section, time_2nd_section, rtime_1st_section, rtime_2nd_section
real(kind=WP)     :: tdiff, time_2nd_gt_1st_section_max, rtime_par_sections

! kh 23.02.21 used for explicit instrumentation and profiling of asynchronous iceberg computations
real(kind=WP), allocatable,dimension(:)     :: time_1st_sections(:), time_2nd_sections(:) !, time_buffer(:)
real(kind=WP), allocatable,dimension(:)     :: time_sections_seq_max(:), time_sections_par_max(:)
real(kind=WP), allocatable,dimension(:, :)  :: time_1st_sections_per_rank(:, :), time_2nd_sections_per_rank(:, :)
real(kind=WP), allocatable,dimension(:)     :: time_saved_by_par(:)
real(kind=WP)     :: time_seq, time_par, total_time_saved_by_par
real(kind=WP)     :: total_time_seq, total_time_par

integer           :: ir, js
integer           :: status(MPI_STATUS_SIZE)

real(kind=real32) :: rtime_setup_mesh, rtime_setup_ocean, rtime_setup_forcing 
real(kind=real32) :: rtime_setup_ice,  rtime_setup_other, rtime_setup_restart

! kh 10.02.21
real(kind=real32) :: mean_rtime(19), max_rtime(19), min_rtime(19)

real(kind=real32) :: runtime_alltimesteps

! kh 02.02.21 
logical             :: first_section_done
logical             :: still_waiting
logical             :: bBreak

real(kind=WP)       :: icb_wait_iterations_counter
real(kind=WP)       :: icb_wait_iterations_counter_min
real(kind=WP)       :: icb_wait_iterations_counter_max
logical             :: ib_async_1st_section_first = .true.
logical             :: ib_async_2nd_section_first = .true.
integer             :: mype_copy
integer             :: limit_list_mype
!integer             :: thread_support_level_required
integer             :: thread_support_level_provided
logical             :: mpi_init_done
integer             :: fesom_thread_info

logical             :: bIcbCalcCycleCompleted

type(t_mesh),             target, save :: mesh

! kh 26.03.21 get current values for ib_async_mode and thread_support_level_required
    call read_namelist_icebergs

#ifndef __oifs
    !ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular dependencies between modules (cpl_driver and g_PARSUP)
    !OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS

    if(ib_async_mode > 0) then
! kh 26.03.21 thread_support_level_required is handled as namelist entry now
!       thread_support_level_required = MPI_THREAD_MULTIPLE
!       thread_support_level_required = MPI_THREAD_SERIALIZED
        call MPI_init_thread(thread_support_level_required, thread_support_level_provided, MPIERR)
!       write(*,*) 'MPI_init_thread: thread_support_level_required, thread_support_level_provided, MPIERR', &
!                   thread_support_level_required, thread_support_level_provided, MPIERR

        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_init_thread(thread_support_level_required, ... failed'
            call par_ex
            stop
        end if

        if (thread_support_level_provided < thread_support_level_required) then
            write(*,*) 'Error: thread_support_level_provided < thread_support_level_required', thread_support_level_provided, thread_support_level_required
            call par_ex
            stop
        end if
    else

! kh 05.08.21 iceberg version
!       call MPI_INIT(i)

! kh 05.08.21 wiso version
        call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, provided, i)
    end if
#endif
    

#if defined (__oasis)
    call cpl_oasis3mct_init(MPI_COMM_FESOM)
#endif
    t1 = MPI_Wtime()

    call par_init

! kh 26.03.21
    if(ib_async_mode > 0) then
        if (mype==0) then
            write(*,*) 'MPI_init_thread: thread_support_level_required, thread_support_level_provided, MPIERR', &
                        thread_support_level_required, thread_support_level_provided, MPIERR
        end if
    end if

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
    !call mesh_setup
    call mesh_setup(mesh)

    if (mype==0) write(*,*) 'FESOM mesh_setup... complete'

! kh 23.02.21 explicit profiling support of asynchronous iceberg computations
    allocate (time_1st_sections(nsteps), time_2nd_sections(nsteps)) !, time_buffer(nsteps))
    allocate (time_sections_seq_max(nsteps), time_sections_par_max(nsteps))
    allocate (time_1st_sections_per_rank(nsteps, npes), time_2nd_sections_per_rank(nsteps, npes))
    allocate (time_saved_by_par(nsteps))

    time_1st_sections           = 0.0_WP
    time_2nd_sections           = 0.0_WP
    time_1st_sections_per_rank  = 0.0_WP
    time_2nd_sections_per_rank  = 0.0_WP
    time_saved_by_par           = 0.0_WP

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
    call forcing_setup(mesh)
    if (mype==0) t4=MPI_Wtime()
    if (use_ice) then 
        call ice_setup(mesh)
        ice_steps_since_upd = ice_ave_steps-1
        ice_update=.true.
        if (mype==0) write(*,*) 'EVP scheme option=', whichEVP
    endif
    if (mype==0) t5=MPI_Wtime()
    call compute_diagnostics(0, mesh) ! allocate arrays for diagnostic
#if defined (__oasis)
    call cpl_oasis3mct_define_unstr(mesh)
    if(mype==0)  write(*,*) 'FESOM ---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv
#endif
   
    !LA: allocate arrays for icebergs
    if (use_icebergs) then
        call allocate_icb
    endif

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
    rtime_fullice               = 0._WP
    rtime_write_restart         = 0._WP
    rtime_write_means           = 0._WP
    rtime_compute_diag          = 0._WP
    rtime_icb_calc              = 0._WP
    rtime_icb_write             = 0._WP

! kh 10.02.21
    time_1st_section            = 0._WP
    time_2nd_section            = 0._WP
    rtime_1st_section           = 0._WP
    rtime_2nd_section           = 0._WP
    time_2nd_gt_1st_section_max = 0._WP
    rtime_par_sections          = 0._WP

    time_1st_sections           = 0._WP
    time_2nd_sections           = 0._WP

! kh 04.02.21 set limit_list_mype for debug purposes
!   limit_list_mype = 9
!   limit_list_mype = 1
    limit_list_mype = -1 ! kh 04.02.21 a value < 0 means "off"

! kh 10.02.21 duplicate communicator to be used in parallel section for async iceberg computations based on OpenMP
    if (ib_async_mode > 0) then
        call MPI_comm_dup(MPI_COMM_FESOM, MPI_COMM_FESOM_IB, MPIERR)
        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_comm_dup(MPI_COMM_FESOM, MPI_COMM_FESOM_IB, MPIERR) failed'
            call par_ex
            stop
        end if

! kh 21.02.21
        call MPI_info_create(fesom_thread_info, MPIERR)
        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_info_create(fesom_thread_info, MPIERR) failed'
            call par_ex
            stop
        end if

        call MPI_info_set(fesom_thread_info, 'thread_id', '0', MPIERR)
        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_info_set(fesom_thread_info, ... failed'
            call par_ex
            stop
        end if
        call MPI_comm_set_info(MPI_COMM_FESOM, fesom_thread_info, MPIERR)
        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_comm_set_info(MPI_COMM_FESOM, ... failed'
            call par_ex
            stop
        end if
        call MPI_info_set(fesom_thread_info, 'thread_id', '1', MPIERR)
        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_info_set(fesom_thread_info, ... failed'
            call par_ex
            stop
        end if
        call MPI_comm_set_info(MPI_COMM_FESOM_IB, fesom_thread_info, MPIERR)
        if (MPIERR /= MPI_SUCCESS) then
            write(*,*) 'MPI_comm_set_info(MPI_COMM_FESOM, ... failed'
            call par_ex
            stop
        end if
    else
        MPI_COMM_FESOM_IB = MPI_COMM_FESOM
    end if
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

! kh 09.02.21
    icb_wait_iterations_counter_min = 1.e20_WP
    icb_wait_iterations_counter_max = 0._WP

! kh 04.02.21 debug output
    if (mype==0) then
        write (*,*) 'ib_async_mode, initial omp_num_threads ', ib_async_mode, omp_get_num_threads()
        write (*,*) 'nsteps, steps_per_ib_step, icb_outfreq :', nsteps, steps_per_ib_step, icb_outfreq
    end if
   
    !___MODEL TIME STEPPING LOOP________________________________________________

    bIcbCalcCycleCompleted = .false.

    n = 1
    do while (n <= nsteps)

        call loop_start_part(n)

        n_ib         = n
        u_wind_ib    = u_wind
        v_wind_ib    = v_wind
        u_ice_ib     = u_ice
        v_ice_ib     = v_ice
        m_ice_ib     = m_ice
        a_ice_ib     = a_ice

! kh 08.03.21 support of different ocean ice and iceberg steps:
! if steps_per_ib_step is configured greater 1 then UV is modified via call oce_timestep_ale(n) -> call update_vel while
! the same asynchronous iceberg computation is still active
        UV_ib        = UV

! kh 15.03.21 support of different ocean ice and iceberg steps:
! if steps_per_ib_step is configured greater 1 then tr_arr is modified via call oce_timestep_ale(n) -> call solve_tracers_ale() while
! the same asynchronous iceberg computation is still active
        tr_arr_ib    = tr_arr

! kh 15.03.21 support of different ocean ice and iceberg steps:
! if steps_per_ib_step is configured greater 1 then Tsurf and Ssurf might be changed while
! the same asynchronous iceberg computation is still active
        Tsurf_ib     = Tsurf
        Ssurf_ib     = Ssurf

! kh 18.03.21 support of different ocean ice and iceberg steps:
! if steps_per_ib_step is configured greater 1 then zbar_3d_n and eta_n might be changed while
! the same asynchronous iceberg computation is still active
        zbar_3d_n_ib = zbar_3d_n
        eta_n_ib     = eta_n

! kh 16.03.21 not modified during overlapping ocean/ice and iceberg computations
!       coriolis_ib      = coriolis
!       coriolis_node_ib = coriolis_node

! kh 02.02.21 check iceberg computations mode:
! ib_async_mode == 0: original sequential behavior for both ice sections (for testing purposes, creating reference results etc.)
! ib_async_mode == 1: OpenMP code active to overlapped computations in first (ocean ice) and second (icebergs) parallel section
! ib_async_mode == 2: OpenMP code active, but computations still serialized via spinlock (for testing purposes)
        if (ib_async_mode == 0) then ! kh 01.03.21 original sequential behavior for both ice sections

! kh 10.03.21 it is not the start of a real parallel section here, but the value of the timer is still of interest
            t1_par_sections = MPI_Wtime()

            call compute_vel_nodes(mesh)

! kh 08.03.21 t1 moved up to here to include the iceberg computation time (like in former FESOM2 paleodyn_icb versions)
           t1 = MPI_Wtime()

! kh 03.03.21 pseudo overlap ocean ice and iceberg calculation to imitate the behavior of the parallel modes (i.e. ib_async_mode > 0)
!           if (use_icebergs .and. mod(n, steps_per_ib_step)==0.0) then
            if (use_icebergs .and. mod(n - 1, steps_per_ib_step)==0) then
                if (mype==0) write(*,*) '*** step n=',n
                t1_icb = MPI_Wtime()
                write(*,*) "LA DEBUG: start iceberg_calculation"
                call iceberg_calculation(mesh,n)
                write(*,*) "LA DEBUG: finish iceberg_calculation"
                t2_icb = MPI_Wtime()
            end if
            !___model sea-ice step__________________________________________________

! kh 08.03.21 t1 moved up to include the iceberg calculation time (like in former FESOM2 paleodyn_icb versions)
!           t1 = MPI_Wtime()
            if(use_ice) then
                call ocean2ice(mesh)
                call update_atm_forcing(n, mesh)
                
                if (ice_steps_since_upd>=ice_ave_steps-1) then
                    ice_update=.true.
                    ice_steps_since_upd = 0
                else
                    ice_update=.false.
                    ice_steps_since_upd=ice_steps_since_upd+1
                endif
                
                if (ice_update) call ice_timestep(n, mesh)
                
                call oce_fluxes_mom(mesh) ! momentum only
                call oce_fluxes(mesh)
            end if  
                
            !###################################
            ! LA check wheather this needs to go inside omp
            call before_oce_step(mesh)
            !###################################

            if (use_icebergs .and. mod(n, steps_per_ib_step)==0.0) then
!               t1_icb = MPI_Wtime()
                !call iceberg_calculation(n)

! kh 08.03.21 add time for call icb2fesom to the end of t2_icb (i.e. time is calculated like in former FESOM2 paleodyn_icb  versions, also see above)
                t1b_icb = MPI_Wtime()
                call icb2fesom(mesh)
!               write(*,*) '*** MASS BALANCE ***'
!               write(*,*) '*** integrated BV: ',SUM(fwbv_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** integrated B: ',SUM(fwb_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** integrated L: ',SUM(fwl_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** integrated E: ',SUM(fwe_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** TOTAL: ',(SUM(fwbv_flux_ib)+SUM(fwb_flux_ib)+SUM(fwl_flux_ib)+SUM(fwe_flux_ib))*dt*steps_per_ib_step
                !alles auf null setzen
                t2b_icb = MPI_Wtime()
!               t2_icb = MPI_Wtime()
                t2_icb = t2_icb + t2b_icb - t1b_icb
                bIcbCalcCycleCompleted = .true.
            end if

! kh 10.03.21 it is not the end of a real parallel section here, but the value of the timer is still of interest
            t2_par_sections = MPI_Wtime()

        else if (ib_async_mode > 0) then ! kh 02.02.21 asynchronous behavior

            t1_par_sections = MPI_Wtime()

            first_section_done = .false.
            mype_copy = mype

! kh 22.02.21
!$omp parallel sections default (none) &
!$omp&  shared(first_section_done, ib_async_1st_section_first, ib_async_2nd_section_first, mype_copy, limit_list_mype, &
!$omp&  mype, &
!$omp&  n, n_ib, mesh, nsteps, steps_per_ib_step, t1, t1_icb, t2_icb, use_ice, use_icebergs, ice_steps_since_upd, ice_ave_steps, ice_update, &
!$omp&  t1_1st_section, t1_2nd_section, t2_1st_section, t2_2nd_section, &
!$omp&  time_1st_section, time_2nd_section, t1_par_sections, t2_par_sections, rtime_par_sections, &
!$omp&  rtime_fullice, rtime_icb_calc, rtime_1st_section, rtime_2nd_section, rtime_compute_diag, rtime_write_means, rtime_write_restart, &
!$omp&  ib_async_mode, icb_wait_iterations_counter_min, icb_wait_iterations_counter_max, thread_support_level_required, bIcbCalcCycleCompleted) &
!$omp&  private(nsub, bBreak, icb_wait_iterations_counter, still_waiting) &
!$omp&  num_threads(2)

!$omp section
! kh 02.02.21 first parallel section

            t1_1st_section = MPI_Wtime()

            if (ib_async_1st_section_first) then
                if (mype_copy <= limit_list_mype) then
                    write(*,*)'mype, mype_copy, thread rank 1st section, omp_num_threads: ', mype, mype_copy, omp_get_thread_num(), omp_get_num_threads()
                end if    
                ib_async_1st_section_first = .false.
            end if

            if (mype_copy <= limit_list_mype) then
                write(*,*) 'now starting ocean ice computatins..., n:', n
            end if

            nsub   = 1
            bBreak = .false.
            do while ((nsub <= steps_per_ib_step) .and. .not. bBreak)
                if (nsub > 1) then
                    call loop_start_part(n)
                end if

! kh 26.03.21 alternatively, a preprocessor definition could be used here
                    if (thread_support_level_required == MPI_THREAD_SERIALIZED) then
                        !$omp critical 
                        call compute_vel_nodes(mesh)
                        !$omp end critical
                    else
                        call compute_vel_nodes(mesh)
                    end if

                !___model sea-ice step__________________________________________________

! kh 08.03.21 do it only once per iceberg step to include the iceberg calculation time (like in former FESOM2 paleodyn_icb versions)
                if (nsub == 1) then
                    t1 = MPI_Wtime()
                end if

                if(use_ice) then

! kh 26.03.21 alternatively, a preprocessor definition could be used here
                    if (thread_support_level_required == MPI_THREAD_SERIALIZED) then
                        !$omp critical 
                        call ocean2ice(mesh)
                        !$omp end critical
                    else
                        call ocean2ice(mesh)
                    end if

! kh 26.03.21 alternatively, a preprocessor defi(mesh)nition could be used here
                    if (thread_support_level_required == MPI_THREAD_SERIALIZED) then
                        !$omp critical 
                        call update_atm_forcing(n, mesh)
                        !$omp end critical
                    else
                        call update_atm_forcing(n, mesh)
                    end if

                    if (ice_steps_since_upd>=ice_ave_steps-1) then
                        ice_update=.true.
                        ice_steps_since_upd = 0
                    else
                        ice_update=.false.
                        ice_steps_since_upd=ice_steps_since_upd+1
                    endif

! kh 26.03.21 alternatively, a preprocessor definition could be used here
                    if (thread_support_level_required == MPI_THREAD_SERIALIZED) then
                        !$omp critical 
                        if (ice_update) call ice_timestep(n, mesh)
                        !$omp end critical
                    else
                        if (ice_update) call ice_timestep(n, mesh)
                    end if
                    
                    call oce_fluxes_mom(mesh) ! momentum only

! kh 26.03.21 alternatively, a preprocessor definition could be used here
                    if (thread_support_level_required == MPI_THREAD_SERIALIZED) then
                        !$omp critical 
                        call oce_fluxes(mesh)
                        !$omp end critical
                    else
                        call oce_fluxes(mesh)
                    end if

                end if ! (use_ice)

                !###################################
                ! LA check wheather this needs to go inside omp
                call before_oce_step(mesh)
                !###################################

                if (nsub < steps_per_ib_step) then
                    if (n == nsteps) then
                        bBreak = .true.
                    else 

! kh 26.03.21 alternatively, a preprocessor definition could be used here
                        if (thread_support_level_required == MPI_THREAD_SERIALIZED) then
                            !$omp critical 
                            call loop_end_part (mesh, .false., bIcbCalcCycleCompleted, n, t1, time_1st_section, time_2nd_section, &
                                                t1_icb, t2_icb, t1_par_sections, t2_par_sections, &
                                                rtime_fullice, rtime_compute_diag, rtime_write_means, rtime_write_restart, rtime_icb_calc, & 
                                                rtime_1st_section, rtime_2nd_section, rtime_par_sections)
                            !$omp end critical
                        else
                            call loop_end_part (mesh, .false., bIcbCalcCycleCompleted, n, t1, time_1st_section, time_2nd_section, &
                                                t1_icb, t2_icb, t1_par_sections, t2_par_sections, &
                                                rtime_fullice, rtime_compute_diag, rtime_write_means, rtime_write_restart, rtime_icb_calc, & 
                                                rtime_1st_section, rtime_2nd_section, rtime_par_sections)
                        end if

                        n = n + 1
                    end if
                end if 

                nsub = nsub + 1
            end do ! ((nsub <= steps_per_ib_step) .and. .not. bBreak)

            if (mype_copy <= limit_list_mype) then
                write(*,*) 'leaving 1st section now, n:', n
            end if

            first_section_done = .true.

            if (ib_async_mode == 2) then
                !$omp flush (first_section_done)
            end if

            t2_1st_section = MPI_Wtime()

!$omp section
! kh 02.02.21 second parallel section

            if (ib_async_2nd_section_first) then
                if (mype_copy <= limit_list_mype) then
                    write(*,*)'mype, mype_copy, thread rank 2nd section, omp_num_threads: ', mype, mype_copy, omp_get_thread_num(), omp_get_num_threads()
                end if
                ib_async_2nd_section_first = .false.
            end if

            if (ib_async_mode == 2) then
! kh 02.02.21 ib_async_mode 2 forces sequential behavior while OpenMP code is already active (e.g. for testing purposes, 
! checking for bit identicality of results (e.g. restart files), ...)

! kh 02.02.21 use busy waiting here for the sake of simplicity (the cores are reserved anyway)
                icb_wait_iterations_counter = 0._WP ! kh 02.02.21 counter to support debugging
                still_waiting = .true. ! kh 02.02.21 only used for debugging

                !$omp flush (first_section_done)
                do while (.not. first_section_done)
                    icb_wait_iterations_counter = icb_wait_iterations_counter + 1.0_WP
                    !$omp flush (first_section_done)
                end do ! (.not. first_section_done)

                still_waiting = .false. ! kh 02.02.21 only used for debugging

                if (mype_copy <= limit_list_mype) then
                    write(*,*) 'icb_wait_iterations_counter: ', icb_wait_iterations_counter
                end if

! kh 05.02.21
                if (icb_wait_iterations_counter < icb_wait_iterations_counter_min) then
                    icb_wait_iterations_counter_min = icb_wait_iterations_counter
                    if (mype_copy <= limit_list_mype) then
                        write (*,*) 'fresh icb_wait_iterations_counter_min, n, n_ib:', icb_wait_iterations_counter_min, n, n_ib
                    end if
                end if

                if (icb_wait_iterations_counter > icb_wait_iterations_counter_max) then
                    icb_wait_iterations_counter_max = icb_wait_iterations_counter
                    if (mype_copy <= limit_list_mype) then
                        write (*,*) 'fresh icb_wait_iterations_counter_max, n, n_ib:', icb_wait_iterations_counter_max, n, n_ib
                    end if
                end if

            end if ! (ib_async_mode == 2) then

! kh 02.02.21 second parallel section (use only net time, i.e. the potential waiting time when using async_mode = 2 is not contained)
            t1_2nd_section = MPI_Wtime()

            if (use_icebergs) then
                if (mype_copy == 0) write(*,*) '*** step n, nib:', n, n_ib

                if (mype_copy <= limit_list_mype) then
                    write(*,*) 'now starting iceberg computations..., n, n_ib:', n, n_ib
                end if
                t1_icb = MPI_Wtime()

! kh 09.03.21
                call iceberg_calculation(mesh,n_ib)
!               call icb2fesom
!               t2_icb = MPI_Wtime()
            end if
!           t2 = MPI_Wtime()

            if (mype_copy <= limit_list_mype) then
                write(*,*) 'leaving 2nd section now, n, n_ib: ', n, n_ib
            end if

            t2_2nd_section = MPI_Wtime()

!$omp end parallel sections

            t2_par_sections = MPI_Wtime()

            time_1st_section = t2_1st_section - t1_1st_section
            time_2nd_section = t2_2nd_section - t1_2nd_section

            time_1st_sections(n_ib) = time_1st_section
            time_2nd_sections(n_ib) = time_2nd_section

            if (mype_copy <= limit_list_mype) then
                write(*,*) 'time_2nd_section, time_1st_section, n, n_ib: ', time_2nd_section, time_1st_section, n, n_ib
            end if
            
            if (time_2nd_section > time_1st_section) then

! kh 09.02.21 suboptimal overlap of 1st section (ocean ice) and 2nd section (icebergs) is indicated

                tdiff = time_2nd_section - time_1st_section

                if (tdiff > time_2nd_gt_1st_section_max) then
                    time_2nd_gt_1st_section_max = tdiff
                    if (mype_copy <= limit_list_mype) then
                        write(*,*) 'time_2nd_section > time_1st_section, new max diff, n, n_ib, t1_1st_section, t2_1st_section, t1_2nd_section, t2_2nd_section: ', &
                                    time_2nd_section, time_1st_section, tdiff, n, n_ib, t1_1st_section, t2_1st_section, t1_2nd_section, t2_2nd_section
                    end if
                else
                    if (mype_copy <= limit_list_mype) then
                        write(*,*) 'time_2nd_section > time_1st_section, diff, n, n_ib: ', time_2nd_section, time_1st_section, tdiff, n, n_ib
                    end if
                endif
            end if

! kh 02.02.21 non overlapping second part of iceberg computations 
            if (use_icebergs) then
!               t1_icb = MPI_Wtime()
!               call iceberg_calculation(n)
                call icb2fesom(mesh)
!               write(*,*) '*** MASS BALANCE ***'
!               write(*,*) '*** integrated BV: ',SUM(fwbv_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** integrated B: ',SUM(fwb_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** integrated L: ',SUM(fwl_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** integrated E: ',SUM(fwe_flux_ib)*dt*steps_per_ib_step
!               write(*,*) '*** TOTAL: ',(SUM(fwbv_flux_ib)+SUM(fwb_flux_ib)+SUM(fwl_flux_ib)+SUM(fwe_flux_ib))*dt*steps_per_ib_step
                !alles auf null setzen
                t2_icb = MPI_Wtime()
                bIcbCalcCycleCompleted = .true.
            end if

! kh 10.02.21 t2 factored out (also see below)
!           t2 = MPI_Wtime()
 
        else
            write(*,*) 'ib_async_mode < 0 is not supported: ', ib_async_mode 
        end if ! (ib_async_mode == 0) then

        call loop_end_part (mesh, .true., bIcbCalcCycleCompleted, n, t1, time_1st_section, time_2nd_section, &
                            t1_icb, t2_icb, t1_par_sections, t2_par_sections, &
                            rtime_fullice, rtime_compute_diag, rtime_write_means, rtime_write_restart, rtime_icb_calc, & 
                            rtime_1st_section, rtime_2nd_section, rtime_par_sections)

        n = n + 1

    end do ! (n <= nsteps)

    if (use_icebergs) then
        t3_icb = MPI_Wtime()

        call iceberg_out
        !call reset_ib_fluxes        
        t4_icb = MPI_Wtime()

        rtime_icb_write = rtime_icb_write + t4_icb - t3_icb
    end if

    if(mype==0) then
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> FESOM STARTS TIME LOOP                                 '//achar(27)//'[0m'
    end if
    !___MODEL TIME STEPPING LOOP________________________________________________
    if (use_global_tides) then
       call foreph_ini(yearnew, month)
    end if

    !do n=1, nsteps        
    !    if (use_global_tides) then
    !       call foreph(mesh)
    !    end if
    !    mstep = n
    !    if (mod(n,logfile_outfreq)==0 .and. mype==0) then
    !        write(*,*) 'FESOM ======================================================='
!   !          write(*,*) 'FESOM step:',n,' day:', n*dt/24./3600.,
    !        write(*,*) 'FESOM step:',n,' day:', daynew,' year:',yearnew 
    !        write(*,*)
    !    end if
!#if !defined (__oifs) || defined (__oasis)
    !        seconds_til_now=INT(dt)*(n-1)
!#end!if
    !    !call clock
    !    
    !    !___compute horizontal velocity on nodes (originaly on elements)________
    !    !call compute_vel_nodes(mesh)
    !    
    !    !___model sea-ice step__________________________________________________
    !    !t1 = MPI_Wtime()
    !    !if(use_ice) then
    !    !    !___compute fluxes from ocean to ice________________________________
    !    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call ocean2ice(n)'//achar(27)//'[0m'
    !    !    !call ocean2ice(mesh)
    !    !    
    !    !    !___compute update of atmospheric forcing____________________________
    !    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing(n)'//achar(27)//'[0m'
    !    !    !t0_frc = MPI_Wtime()
    !    !    !call update_atm_forcing(n, mesh)
    !    !    !t1_frc = MPI_Wtime()            
    !    !    !___compute ice step________________________________________________
    !    !    !if (ice_steps_since_upd>=ice_ave_steps-1) then
    !    !    !    ice_update=.true.
    !    !    !    ice_steps_since_upd = 0
    !    !    !else
    !    !    !    ice_update=.false.
    !    !    !    ice_steps_since_upd=ice_steps_since_upd+1
    !    !    !endif
    !    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
    !    !    !if (ice_update) call ice_timestep(n, mesh)  
    !    !    !!___compute fluxes to the ocean: heat, freshwater, momentum_________
    !    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call oce_fluxes_mom...'//achar(27)//'[0m'
    !    !    !call oce_fluxes_mom(mesh) ! momentum only
    !    !    !call oce_fluxes(mesh)
    !    !end if
    !    !call before_oce_step(mesh) ! prepare the things if required
    !    !t2 = MPI_Wtime()
    !    
    !    !!___model ocean step____________________________________________________
    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call oce_timestep_ale'//achar(27)//'[0m'
    !    !call oce_timestep_ale(n, mesh)
    !    !t3 = MPI_Wtime()
    !    !!___compute energy diagnostics..._______________________________________
    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics(1, mesh)'//achar(27)//'[0m'
    !    !call compute_diagnostics(1, mesh)
    !    !t4 = MPI_Wtime()
    !    !___prepare output______________________________________________________
    !    !if (flag_debug .and. mype==0)  print *, achar(27)//'[34m'//' --> call output (n)'//achar(27)//'[0m'
    !    !call output (n, mesh)
    !    !t5 = MPI_Wtime()
    !    !call restart(n, .false., .false., mesh)
    !    !t6 = MPI_Wtime()
    !    !
    !    !rtime_fullice       = rtime_fullice       + t2 - t1
    !    !rtime_compute_diag  = rtime_compute_diag  + t4 - t3
    !    !rtime_write_means   = rtime_write_means   + t5 - t4   
    !    !rtime_write_restart = rtime_write_restart + t6 - t5
    !    !rtime_read_forcing  = rtime_read_forcing  + t1_frc - t0_frc
    !end do
    
    call finalize_output()
    
    !___FINISH MODEL RUN________________________________________________________

    call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
    if (mype==0) then
       t1 = MPI_Wtime()
       runtime_alltimesteps = real(t1-t0,real32)
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
    mean_rtime(14) = rtime_icb_calc
    mean_rtime(15) = rtime_icb_write

! kh 10.02.21
    mean_rtime(16) = rtime_1st_section
    mean_rtime(17) = rtime_2nd_section
    mean_rtime(18) = rtime_par_sections

    max_rtime(1:18) = mean_rtime(1:18)
    min_rtime(1:18) = mean_rtime(1:18)

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 18, MPI_REAL, MPI_SUM, MPI_COMM_FESOM, MPIerr)
    mean_rtime(1:18) = mean_rtime(1:18) / real(npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  18, MPI_REAL, MPI_MAX, MPI_COMM_FESOM, MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  18, MPI_REAL, MPI_MIN, MPI_COMM_FESOM, MPIerr)

! kh 23.02.21 start collecting profiling data of omp sections
    if(mype == 0) then
        time_1st_sections_per_rank(:, 1) = time_1st_sections
        time_2nd_sections_per_rank(:, 1) = time_2nd_sections
        do ir = 1, npes - 1
            call MPI_recv(time_1st_sections_per_rank(:, ir + 1), nsteps, MPI_DOUBLE, ir, 0, MPI_COMM_FESOM, status, MPIerr)
        end do
        do ir = 1, npes - 1
            call MPI_recv(time_2nd_sections_per_rank(:, ir + 1), nsteps, MPI_DOUBLE, ir, 0, MPI_COMM_FESOM, status, MPIerr)
        end do
    else
        call MPI_send(time_1st_sections, nsteps, MPI_DOUBLE, 0, 0, MPI_COMM_FESOM, MPIerr)
        call MPI_send(time_2nd_sections, nsteps, MPI_DOUBLE, 0, 0, MPI_COMM_FESOM, MPIerr)
    end if

    if (ib_async_mode > 0) then
        if(mype == 0) then
            do js = 1, nsteps
                time_sections_seq_max(js) = 0.0_WP
                time_sections_par_max(js) = 0.0_WP
                do ir = 0, npes - 1
                    time_seq = time_1st_sections_per_rank(js, ir + 1) + time_2nd_sections_per_rank(js, ir + 1)
                    if (time_1st_sections_per_rank(js, ir + 1) > time_2nd_sections_per_rank(js, ir + 1)) then
                        time_par = time_1st_sections_per_rank(js, ir + 1)
                    else
                        time_par = time_2nd_sections_per_rank(js, ir + 1)
                    end if
                    if(time_seq > time_sections_seq_max(js)) then
                        time_sections_seq_max(js) = time_seq
                    end if
                    if(time_par > time_sections_par_max(js)) then
                        time_sections_par_max(js) = time_par
                    end if
                end do
            end do

            total_time_seq = 0.0_WP
            total_time_par = 0.0_WP

            total_time_saved_by_par = 0.0_WP
            do js = 1, nsteps
                time_saved_by_par(js) = time_sections_seq_max(js) - time_sections_par_max(js)

! kh 24.02.21 assert time_saved_by_par(js) >= 0.0
                if(time_saved_by_par(js) < 0.0) then
                    write(*,*) 'Error: time_saved_by_par(js) < 0.0', time_saved_by_par
                end if

                total_time_saved_by_par = total_time_saved_by_par + time_saved_by_par(js)

                total_time_seq = total_time_seq + time_sections_seq_max(js)
                total_time_par = total_time_par + time_sections_par_max(js)

! kh 23.03.21 write only the last line
                if (js == nsteps) then
                    write(*,*) 'time_saved_by_par(js), time_sections_seq_max(js), time_sections_par_max(js), js, total_time_saved_by_par:', &
                               time_saved_by_par(js), time_sections_seq_max(js), time_sections_par_max(js), js, total_time_saved_by_par
                end if
            end do

            write(*,*) 'total_time_seq, total_time_par, total_time_saved_by_par', total_time_seq, total_time_par, total_time_saved_by_par
            write(*,*)
        end if ! (mype == 0) then
    end if !(ib_async_mode > 0) then
! kh 23.02.21 end collecting profiling data of omp sections
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
        write(*,*) '  runtime icb calc:', mean_rtime(14), min_rtime(14), max_rtime(14)
        write(*,*) '  runtime icb write:', mean_rtime(15), min_rtime(15), max_rtime(15)

        if (ib_async_mode > 0) then
            write(*,*) '  runtime par. sections:', mean_rtime(18), min_rtime(18), max_rtime(18)

            write(*,*) '  runtime 1st par. section:', mean_rtime(16), min_rtime(16), max_rtime(16)
            write(*,*) '  runtime 2nd par. section:', mean_rtime(17), min_rtime(17), max_rtime(17)
        else
            write(*,*) '  runtime potential par. sections:', mean_rtime(18), min_rtime(18), max_rtime(18)
        end if
        
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

! kh 19.03.21
    deallocate (time_1st_sections, time_2nd_sections) !, time_buffer)
    deallocate (time_sections_seq_max, time_sections_par_max)
    deallocate (time_1st_sections_per_rank, time_2nd_sections_per_rank)
    deallocate (time_saved_by_par)

! kh 10.02.21 free duplicated communicator used in parallel section for async iceberg computations based on OpenMP
! kh 11.02.21 not meaningful here
!   if (ib_async_mode > 0) then
!       call MPI_Comm_free (MPI_COMM_FESOM_IB, MPIERR)
!   end if

    call par_ex
end program main



! kh 25.02.21 start part from main loop factored out for the async iceberg implementation
subroutine loop_start_part(n)
    use o_PARAM     ! mstep
    use g_config    ! logfile_outfreq, dt
    use g_PARSUP    ! mype
    use g_clock     ! daynew, yearnew
    use cpl_driver  ! seconds_til_now

    implicit none
    integer, intent(in)     :: n

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

end subroutine loop_start_part



! kh 25.02.21 end part from main loop factored out for the async iceberg implementation
subroutine loop_end_part (mesh, bOuterLoopCall, bIcbCalcCycleCompleted, n, t1, time_1st_section, time_2nd_section, &
                          t1_icb, t2_icb, t1_par_sections, t2_par_sections, &
                          rtime_fullice, rtime_compute_diag, rtime_write_means, rtime_write_restart, rtime_icb_calc, & 
                          rtime_1st_section, rtime_2nd_section, rtime_par_sections)
    use o_PARAM         ! WP
    use g_PARSUP        ! MPI_Wtime()
    use io_RESTART      ! restart(...)
    use io_MEANDATA     ! output(...)
    use diagnostics     ! compute_diagnostics(...)
    use g_config        ! use_icebergs
    
    use g_parsup
    use MOD_MESH

    implicit none
    logical, intent(in)             :: bOuterLoopCall
    logical, intent(inout)          :: bIcbCalcCycleCompleted
    integer, intent(in)             :: n
    real(kind=WP), intent(in)       :: t1, time_1st_section, time_2nd_section
    real(kind=WP), intent(in)       :: t1_icb, t2_icb, t1_par_sections, t2_par_sections
    real(kind=WP), intent(inout)    :: rtime_fullice, rtime_compute_diag, rtime_write_means, rtime_write_restart, rtime_icb_calc
    real(kind=WP), intent(inout)    :: rtime_1st_section, rtime_2nd_section, rtime_par_sections

    real(kind=WP)                   :: t2, t3, t4, t5, t6
type(t_mesh), intent(in) , target :: mesh
#include "associate_mesh.h"

    t2 = MPI_Wtime()
        
!___model ocean step____________________________________________________
    call oce_timestep_ale(n, mesh)
    t3 = MPI_Wtime()
    call compute_diagnostics(1, mesh)
    t4 = MPI_Wtime()
    !___prepare output______________________________________________________
    call output (n, mesh)
    if (use_icebergs .and. bIcbCalcCycleCompleted) then
        call reset_ib_fluxes
    end if

    t5 = MPI_Wtime()
    call restart(n, .false., .false., mesh)
    t6 = MPI_Wtime()

! kh 08.03.21 do it only once per iceberg step if ib_async_mode > 0 to include the iceberg calculation time (like in former FESOM2 paleodyn_icb versions)
    if (bOuterLoopCall) then
        rtime_fullice = rtime_fullice + t2 - t1
    end if 

    rtime_compute_diag  = rtime_compute_diag  + t4 - t3
    rtime_write_means   = rtime_write_means   + t5 - t4   
    rtime_write_restart = rtime_write_restart + t6 - t5

! kh 08.03.21 do it only once per iceberg step to include the iceberg calculation time (like in former FESOM2 paleodyn_icb versions)
    if (bIcbCalcCycleCompleted) then
        rtime_icb_calc = rtime_icb_calc + t2_icb - t1_icb
    end if

! kh 10.02.21 do it only once per iceberg step if ib_async_mode > 0
    if (bOuterLoopCall) then
        rtime_1st_section = rtime_1st_section + time_1st_section
        rtime_2nd_section = rtime_2nd_section + time_2nd_section

! kh 24.02.21        
        rtime_par_sections = rtime_par_sections + t2_par_sections - t1_par_sections
    end if

    bIcbCalcCycleCompleted = .false.

end subroutine loop_end_part
