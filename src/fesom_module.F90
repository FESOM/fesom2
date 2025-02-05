! synopsis: save any derived types we initialize
!           so they can be reused after fesom_init
module fesom_main_storage_module
  use iceberg_step
  USE MOD_MESH
  USE MOD_ICE
  USE MOD_TRACER
  USE MOD_PARTIT
  USE MOD_PARSUP
  USE MOD_DYN
  USE o_ARRAYS
  USE o_PARAM
  use g_clock
  use g_config
  use g_comm_auto
  use g_forcing_arrays
  use io_RESTART
  use io_MEANDATA
  use io_mesh_info
  use diagnostics
  use mo_tidal
  use tracer_init_interface
  use ocean_setup_interface
  use ice_setup_interface
  use ocean2ice_interface
  use oce_fluxes_interface
  use update_atm_forcing_interface
  use before_oce_step_interface
  use oce_timestep_ale_interface
  use read_mesh_interface
  use fesom_version_info_module
  use command_line_options_module
  use, intrinsic :: iso_fortran_env, only : real32
  use g_forcing_param, only: use_landice_water, use_age_tracer
  use landice_water_init_interface
  use age_tracer_init_interface
  use iceberg_params
  use iceberg_step
  ! Define icepack module

#if defined (__icepack)
  use icedrv_main,          only: set_icepack, init_icepack, alloc_icepack
#endif

#if defined (__oasis)
  use cpl_driver
#endif

! define recom module
#if defined (__recom)
  use recom_init_interface
  use recom_interface
#endif

  implicit none
    
  type :: fesom_main_storage_type

    integer           :: n, from_nstep, offset, row, i, provided
    integer           :: which_readr ! read which restart files (0=netcdf, 1=core dump,2=dtype)
    integer           :: total_nsteps
    integer, pointer  :: mype, npes, MPIerr, MPI_COMM_FESOM, MPI_COMM_WORLD, MPI_COMM_FESOM_IB
#if defined(__recom) && defined(__usetp)
    integer, pointer  :: my_fesom_group, MPI_COMM_FESOM_WORLD, MPI_COMM_FESOM_SAME_RANK_IN_GROUPS
#endif
    real(kind=WP)     :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t0_ice, t1_ice, t0_frc, t1_frc
    real(kind=WP)     :: rtime_fullice,    rtime_write_restart, rtime_write_means, rtime_compute_diag, rtime_read_forcing
    real(kind=real32) :: rtime_setup_mesh, rtime_setup_ocean, rtime_setup_forcing 
    real(kind=real32) :: rtime_setup_ice,  rtime_setup_other, rtime_setup_restart
    real(kind=real32) :: runtime_alltimesteps
#if defined (__recom)
    real(kind=WP)     :: t0_recom, t1_recom
    real(kind=real32) :: rtime_setup_recom, rtime_compute_recom
#endif

    type(t_mesh)   mesh
    type(t_tracer) tracers
    type(t_dyn)    dynamics
    type(t_partit) partit
    type(t_ice)    ice


    character(LEN=256)               :: dump_dir, dump_filename
    logical                          :: L_EXISTS
    type(t_mesh)   mesh_copy
    type(t_tracer) tracers_copy
    type(t_dyn)    dynamics_copy
    type(t_ice)    ice_copy

    character(LEN=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_version_txt
    integer mpi_version_len
    logical fesom_did_mpi_init
    
  end type
  type(fesom_main_storage_type), save, target :: f

end module


! synopsis: main FESOM program split into 3 parts
!           this way FESOM can e.g. be used as a library with an external time loop driver
!           used with IFS-FESOM
module fesom_module
#if defined  __ifsinterface
  use, intrinsic :: ieee_exceptions
#endif
  implicit none
  public fesom_init, fesom_runloop, fesom_finalize
  private

contains
 
  subroutine fesom_init(fesom_total_nsteps)
      use fesom_main_storage_module
#if defined(__MULTIO)
      use iom
#endif
      use cpl_driver
      integer, intent(out) :: fesom_total_nsteps
      ! EO parameters
      logical mpi_is_initialized
      integer              :: tr_num

#if defined(__recom) && defined(__usetp)
! kh 11.11.21 multi FESOM group loop parallelization
      integer             :: npes_fesom_world
      integer             :: mype_fesom_world
      integer             :: processes_per_group
      integer             :: npes_check
      integer             :: mype_check

! kh 26.11.21 get current value for num_fesom_groups
    call read_namelist_run_config
#endif

#if !defined  __ifsinterface
      if(command_argument_count() > 0) then
        call command_line_options%parse()
        stop
      end if
#endif

!SUVI: disable overflow, underflow for entire model when used  in coupled with ifs
!      bad practice, use it only in case of Emergency
!#if defined  __ifsinterface
!    call ieee_set_halting_mode(ieee_overflow, .false.)
!    call ieee_set_halting_mode(ieee_underflow, .false.)
!    call ieee_set_halting_mode(ieee_invalid, .false.)
!#endif

      
      mpi_is_initialized = .false.
      f%fesom_did_mpi_init = .false.

#ifndef __oifs
        !ECHAM6-FESOM2 coupling: cpl_oasis3mct_init is called here in order to avoid circular dependencies between modules (cpl_driver and g_PARSUP)
        !OIFS-FESOM2 coupling: does not require MPI_INIT here as this is done by OASIS
        call MPI_Initialized(mpi_is_initialized, f%i)
        if(.not. mpi_is_initialized) then
            ! TODO: do not initialize MPI here if it has been initialized already, e.g. via IFS when fesom is called as library (__ifsinterface is defined)
            call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, f%provided, f%i)
            f%fesom_did_mpi_init = .true.
        end if
#endif

#if defined (__oasis)
!        call cpl_oasis3mct_init(f%partit,f%partit%MPI_COMM_FESOM)
! kh 02.12.21
#if defined(__recom) && defined(__usetp)
        call cpl_oasis3mct_init(f%partit, f%partit%MPI_COMM_FESOM, num_fesom_groups)
#else
        call cpl_oasis3mct_init(f%partit, f%partit%MPI_COMM_FESOM)
#endif
#endif
        f%t1 = MPI_Wtime()

        call par_init(f%partit)

        f%mype          =>f%partit%mype
        f%MPIerr        =>f%partit%MPIerr
        f%MPI_COMM_FESOM=>f%partit%MPI_COMM_FESOM
        f%MPI_COMM_FESOM_IB=>f%partit%MPI_COMM_FESOM_IB
        f%MPI_COMM_WORLD=>f%partit%MPI_COMM_WORLD

        f%npes          =>f%partit%npes

#if defined(__recom) && defined(__usetp)
! kh 26.11.21 prepare communicator splitting for multi FESOM group loop parallelization
    f%my_fesom_group=>f%partit%my_fesom_group

    f%MPI_COMM_FESOM_WORLD=> f%partit%MPI_COMM_FESOM_WORLD
    f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS=> f%partit%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS

    f%MPI_COMM_FESOM_WORLD = f%MPI_COMM_FESOM
    npes_fesom_world     = f%npes
    mype_fesom_world     = f%mype
    if(mype_fesom_world == 0) then
        write(*,*) 'npes_fesom_world, num_fesom_groups', npes_fesom_world, num_fesom_groups
    end if
    if(mod(npes_fesom_world, num_fesom_groups) /= 0) then
        if(mype_fesom_world == 0) then
            write(*,*) 'MPI_comm_split mismatch npes_fesom_world, num_fesom_groups', npes_fesom_world, num_fesom_groups
        end if
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    processes_per_group = npes_fesom_world / num_fesom_groups
    if(mype_fesom_world == 0) then
        write(*,*) 'processes_per_group', processes_per_group
    end if
    f%npes           = processes_per_group
    f%my_fesom_group = mype_fesom_world / processes_per_group
    f%mype           = mod(mype_fesom_world, processes_per_group)

! kh 26.11.21 split to num_fesom_groups
    call MPI_comm_split(f%MPI_COMM_FESOM_WORLD, f%my_fesom_group, 0, f%MPI_COMM_FESOM, f%MPIerr)
    if (f%MPIerr /= MPI_SUCCESS) then
        write(*,*) 'MPI_comm_split(MPI_COMM_FESOM_WORLD, my_fesom_group, 0, MPI_COMM_FESOM, MPIERR) failed'
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    call MPI_comm_size(f%MPI_COMM_FESOM, npes_check, f%MPIerr)
    if(f%MPIerr /= MPI_SUCCESS) then
        write(*,*) 'MPI_comm_size(MPI_COMM_FESOM, npes_check, MPIERR) failed'
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    call MPI_comm_rank(f%MPI_COMM_FESOM, mype_check, f%MPIerr)
    if(f%MPIerr /= MPI_SUCCESS) then
        write(*,*) 'MPI_comm_rank(MPI_COMM_FESOM, mype_check, MPIERR) failed'
        call par_ex(f%partit%MPI_COMM_FESOM, f%mype)
        stop
    end if

    if(npes_check /= f%npes) then
        write(*,*) 'npes mismatch, npes, npes_check', f%npes, npes_check
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    if(mype_check /= f%mype) then
        write(*,*) 'mype mismatch, mype, mype_check', f%mype, mype_check
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

! kh 17.11.21 group same ranks in each group for broadcasting
!   write(*,*) 'mype, my_fesom_group', mype, my_fesom_group

    call MPI_comm_split(f%MPI_COMM_FESOM_WORLD, f%mype, f%my_fesom_group, f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, f%MPIERR)
    if (f%MPIERR /= MPI_SUCCESS) then
        write(*,*) 'MPI_comm_split(MPI_COMM_FESOM_WORLD, mype, my_fesom_group, MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, MPIERR) failed'
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    call MPI_comm_size(f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, npes_check, f%MPIERR)
    if(f%MPIERR /= MPI_SUCCESS) then
        write(*,*) 'MPI_comm_size(MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, npes_check, MPIERR) failed'
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    call MPI_comm_rank(f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, mype_check, f%MPIERR)
    if(f%MPIERR /= MPI_SUCCESS) then
        write(*,*) 'MPI_comm_rank(MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, mype_check, MPIERR) failed'
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    if(npes_check /= num_fesom_groups) then
        write(*,*) 'npes mismatch, num_fesom_groups, npes_check', num_fesom_groups, npes_check
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

    if(mype_check /= f%my_fesom_group) then
        write(*,*) 'mype mismatch, my_fesom_group, mype_check', f%my_fesom_group, mype_check
        call par_ex(f%MPI_COMM_FESOM, f%mype)
        stop
    end if

! kh 29.02.22
    if(f%my_fesom_group==0) then
#endif

        if(f%mype==0) then
            write(*,*)
            print *,"FESOM2 git SHA: "//fesom_git_sha()
            call MPI_Get_library_version(f%mpi_version_txt, f%mpi_version_len, f%MPIERR)
            print *,"MPI library version: "//trim(f%mpi_version_txt)
            print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
            print *, achar(27)//'[7;32m'//' --> FESOM BUILDS UP MODEL CONFIGURATION                    '//achar(27)//'[0m'
        end if

#if defined(__recom) && defined(__usetp)
    end if
#endif

        !=====================
        ! Read configuration data,  
        ! load the mesh and fill in 
        ! auxiliary mesh arrays
        !=====================
        call setup_model(f%partit)  ! Read Namelists, always before clock_init
        call clock_init(f%partit)   ! read the clock file 
        call get_run_steps(fesom_total_nsteps, f%partit)
        f%total_nsteps=fesom_total_nsteps
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call mesh_setup'//achar(27)//'[0m'
        call mesh_setup(f%partit, f%mesh)

! kh 29.02.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

        if (f%mype==0) write(*,*) 'FESOM mesh_setup... complete'

#if defined(__recom) && defined(__usetp)
        end if
#endif 

        !=====================
        ! Allocate field variables 
        ! and additional arrays needed for 
        ! fancy advection etc.  
        !=====================
#if defined (__oasis)
        !---wiso-code
        IF (lwiso) THEN
          nsend = nsend + 6       ! add number of water isotope tracers to coupling parameter nsend, nrecv
          nrecv = nrecv + 6
        END IF
        !---wiso-code-end
        IF (use_icebergs) THEN
          nrecv = nrecv + 2
        END IF
#endif

        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call check_mesh_consistency'//achar(27)//'[0m'
        call check_mesh_consistency(f%partit, f%mesh)
        if (f%mype==0) f%t2=MPI_Wtime()

        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call dynamics_init'//achar(27)//'[0m'
        call dynamics_init(f%dynamics, f%partit, f%mesh)
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call tracer_init'//achar(27)//'[0m'
        call tracer_init(f%tracers, f%partit, f%mesh)                ! allocate array of ocean tracers (derived type "t_tracer")

        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call arrays_init'//achar(27)//'[0m'
        call arrays_init(f%tracers%num_tracers, f%partit, f%mesh)    ! allocate other arrays (to be refactured same as tracers in the future)
        
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ocean_setup'//achar(27)//'[0m'
        call ocean_setup(f%dynamics, f%tracers, f%partit, f%mesh)

        ! recom setup
#if defined (__recom)
#if defined (__usetp)
        if(f%my_fesom_group==0) then
#endif 
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call recom_init'//achar(27)//'[0m'
#if defined (__usetp)
        end if
#endif

        f%t0_recom=MPI_Wtime()
        call recom_init(f%tracers, f%partit, f%mesh) ! adjust values for recom tracers (derived type "t_tracer")
        f%t1_recom=MPI_Wtime()

#if defined (__usetp)
        if(f%my_fesom_group==0) then
#endif

        if (f%mype==0) write(*,*) 'RECOM recom_init... complete'

#if defined (__usetp)
        end if
#endif

#endif

        if (f%mype==0) then

! kh 29.02.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif

           write(*,*) 'FESOM ocean_setup... complete'

#if defined(__recom) && defined(__usetp)
        end if
#endif

           f%t3=MPI_Wtime()
        endif
        call forcing_setup(f%partit, f%mesh)

        if (f%mype==0) f%t4=MPI_Wtime()
        if (use_ice) then 
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_setup'//achar(27)//'[0m'
            call ice_setup(f%ice, f%tracers, f%partit, f%mesh)
            f%ice%ice_steps_since_upd = f%ice%ice_ave_steps-1
            f%ice%ice_update=.true.

#if defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 
            
            if (f%mype==0) write(*,*) 'EVP scheme option=', f%ice%whichEVP

#if defined(__usetp)
        end if
#endif 

        else 
            ! create a dummy ice derived type with only a_ice, m_ice, m_snow and 
            ! uvice since oce_timesteps still needs in moment
            ! ice as an input for mo_convect(ice, partit, mesh), call 
            ! compute_vel_rhs(ice, dynamics, partit, mesh),  
            ! call write_step_info(...) and call check_blowup(...)
            call ice_init_toyocean_dummy(f%ice, f%partit, f%mesh)
        endif
        
        if (f%mype==0) f%t5=MPI_Wtime()

        call compute_diagnostics(0, f%dynamics, f%tracers, f%ice, f%partit, f%mesh) ! allocate arrays for diagnostic
        
        !---fwf-code-begin
        if(f%mype==0)  write(*,*) 'use_landice_water', use_landice_water
        if(use_landice_water) call landice_water_init(f%partit, f%mesh)
        !---fwf-code-end

        !---age-code-begin
        if(f%mype==0)  write(*,*) 'use_age_tracer', use_age_tracer
        if(use_age_tracer) call age_tracer_init(f%partit, f%mesh)
        !---age-code-end
#if defined (__oasis)

! kh 30.11.21 only mype == 0 in my_fesom_group == 0 handles coupling with extern models
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

        call cpl_oasis3mct_define_unstr(f%partit, f%mesh)

        if(f%mype==0)  write(*,*) 'FESOM ---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv

#if defined(__recom) && defined(__usetp)
        end if
#endif

#if defined(__recom) && defined(__usetp)
! kh 03.12.21
!    call MPI_Barrier(f%MPI_COMM_FESOM_WORLD, f%MPIERR)
    if(num_fesom_groups > 1) then

        call MPI_Bcast(cpl_send, sizeof(cpl_send), MPI_CHARACTER, 0, f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, f%MPIerr)
        call MPI_Bcast(cpl_recv, sizeof(cpl_recv), MPI_CHARACTER, 0, f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, f%MPIerr)

! kh 10.12.21 needed in SUBROUTINE net_rec_from_atm(action)
        call MPI_Bcast(target_root,             1, MPI_INTEGER,   0, f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, f%MPIerr)
    end if
#endif

#endif
    
        ! --------------
        ! LA icebergs: 2023-05-17 
        if (use_icebergs) then
            call allocate_icb(f%partit)
        endif
        ! --------------

#if defined (__icepack)
        !=====================
        ! Setup icepack
        !=====================
        if (f%mype==0) write(*,*) 'Icepack: reading namelists from namelist.icepack'
        call set_icepack(f%ice, f%partit)
        call alloc_icepack
        call init_icepack(f%ice, f%tracers%data(1), f%mesh)
        if (f%mype==0) write(*,*) 'Icepack: setup complete'
#endif
        call clock_newyear                        ! check if it is a new year
        if (f%mype==0) f%t6=MPI_Wtime()
        !___CREATE NEW RESTART FILE IF APPLICABLE___________________________________
        call restart(0, 0, 0, r_restart, f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        if (f%mype==0) f%t7=MPI_Wtime()
        ! store grid information into netcdf file

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif

        if (.not. r_restart) call write_mesh_info(f%partit, f%mesh)

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif

        !___IF RESTART WITH ZLEVEL OR ZSTAR IS DONE, ALSO THE ACTUAL LEVELS AND ____
        !___MIDDEPTH LEVELS NEEDS TO BE CALCULATET AT RESTART_______________________
        if (r_restart .and. .not. f%which_readr==2) then
            call restart_thickness_ale(f%partit, f%mesh)
        end if
        if (f%mype==0) then
           f%t8=MPI_Wtime()
    
           f%rtime_setup_mesh    = real( f%t2 - f%t1              ,real32)
           f%rtime_setup_ocean   = real( f%t3 - f%t2              ,real32)
           f%rtime_setup_forcing = real( f%t4 - f%t3              ,real32)
           f%rtime_setup_ice     = real( f%t5 - f%t4              ,real32)
           f%rtime_setup_restart = real( f%t7 - f%t6              ,real32)
           f%rtime_setup_other   = real((f%t8 - f%t7) + (f%t6 - f%t5) ,real32)

#if defined (__recom)
           f%rtime_setup_recom   = real( f%t1_recom - f%t0_recom  ,real32)
#endif

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif

           write(*,*) '=========================================='
           write(*,*) 'MODEL SETUP took on mype=0 [seconds]      '
           write(*,*) 'runtime setup total      ',real(f%t8-f%t1,real32)      
           write(*,*) ' > runtime setup mesh    ',f%rtime_setup_mesh   
           write(*,*) ' > runtime setup ocean   ',f%rtime_setup_ocean  
           write(*,*) ' > runtime setup forcing ',f%rtime_setup_forcing
           write(*,*) ' > runtime setup ice     ',f%rtime_setup_ice    
           write(*,*) ' > runtime setup restart ',f%rtime_setup_restart
           write(*,*) ' > runtime setup other   ',f%rtime_setup_other
#if defined (__recom)
           write(*,*) ' > runtime setup recom   ',f%rtime_setup_recom
#endif
            write(*,*) '============================================' 

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif

        endif

#if defined(__MULTIO)
          call iom_send_fesom_domains(f%partit, f%mesh)
#endif

    !    f%dump_dir='DUMP/'
    !    INQUIRE(file=trim(f%dump_dir), EXIST=f%L_EXISTS)
    !    if (.not. f%L_EXISTS) call system('mkdir '//trim(f%dump_dir))

    !    write (f%dump_filename, "(A7,I7.7)") "t_mesh.", f%mype
    !    open  (f%mype+300, file=TRIM(f%dump_dir)//trim(f%dump_filename), status='replace', form="unformatted")
    !    write (f%mype+300) f%mesh
    !    close (f%mype+300)

    !    open  (f%mype+300, file=trim(f%dump_filename), status='old', form="unformatted")
    !    read  (f%mype+300) f%mesh_copy
    !    close (f%mype+300)
         
    !    write (f%dump_filename, "(A9,I7.7)") "t_tracer.", f%mype
    !    open  (f%mype+300, file=TRIM(f%dump_dir)//trim(f%dump_filename), status='replace', form="unformatted")
    !    write (f%mype+300) f%tracers
    !    close (f%mype+300)

    !    open  (f%mype+300, file=trim(f%dump_filename), status='old', form="unformatted")
    !    read  (f%mype+300) f%dynamics_copy
    !    close (f%mype+300)

    !    write (f%dump_filename, "(A9,I7.7)") "t_dynamics.", f%mype
    !    open  (f%mype+300, file=TRIM(f%dump_dir)//trim(f%dump_filename), status='replace', form="unformatted")
    !    write (f%mype+300) f%dynamics
    !    close (f%mype+300)

    !    open  (f%mype+300, file=trim(f%dump_filename), status='old', form="unformatted")
    !    read  (f%mype+300) f%tracers_copy
    !    close (f%mype+300)
    
    !call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype)
    !stop
    !         
    !    if (f%mype==10) write(,) f%mesh1%ssh_stiff%values-f%mesh%ssh_stiff%value    
  
    ! Initialize timers
    f%rtime_fullice       = 0._WP
    f%rtime_write_restart = 0._WP
    f%rtime_write_means   = 0._WP
    f%rtime_compute_diag  = 0._WP
    f%rtime_read_forcing  = 0._WP
#if defined (__recom)
    f%rtime_compute_recom = 0._WP
#endif

    f%from_nstep = 1

    !enter mesh and partit data. 
    !$ACC ENTER DATA COPYIN (f) 
    !$ACC ENTER DATA COPYIN (f%mesh, f%mesh%coriolis_node, f%mesh%nn_num, f%mesh%nn_pos) 
    !$ACC ENTER DATA COPYIN (f%mesh%ssh_stiff, f%mesh%ssh_stiff%rowptr) 
    !$ACC ENTER DATA COPYIN (f%mesh%gradient_sca, f%mesh%metric_factor, f%mesh%elem_area, f%mesh%area, f%mesh%edge2D_in) 
    !$ACC ENTER DATA COPYIN (f%mesh%elem2D_nodes, f%mesh%ulevels, f%mesh%ulevels_nod2d, f%mesh%edges, f%mesh%edge_tri) 
    !$ACC ENTER DATA COPYIN (f%partit, f%partit%eDim_nod2D, f%partit%myDim_edge2D) 
    !$ACC ENTER DATA COPYIN (f%partit%myDim_elem2D, f%partit%myDim_nod2D, f%partit%myList_edge2D) 

    !$ACC ENTER DATA COPYIN (f%mesh%elem_cos, f%mesh%edge_cross_dxdy, f%mesh%elem2d_nodes, f%mesh%nl) 
    !$ACC ENTER DATA COPYIN (f%mesh%nlevels_nod2D, f%mesh%nod_in_elem2D, f%mesh%nod_in_elem2D_num) 
    !$ACC ENTER DATA COPYIN (f%mesh%edge_dxdy, f%mesh%nlevels, f%mesh%ulevels_nod2D_max) 
    !$ACC ENTER DATA COPYIN (f%mesh%areasvol, f%mesh%nlevels_nod2D_min) 
    !$ACC ENTER DATA CREATE (f%mesh%helem, f%mesh%hnode, f%mesh%hnode_new, f%mesh%zbar_3d_n, f%mesh%z_3d_n)
    !do n=f%from_nstep, f%from_nstep-1+current_nsteps
    !$ACC ENTER DATA COPYIN  (f%ice)
    !$ACC ENTER DATA CREATE  (f%ice%data, f%ice%work, f%ice%work%fct_massmatrix) 
    !$ACC ENTER DATA CREATE  (f%ice%delta_min, f%ice%Tevp_inv, f%ice%cd_oce_ice) 
    !$ACC ENTER DATA CREATE  (f%ice%work%fct_tmax, f%ice%work%fct_tmin) 
    !$ACC ENTER DATA CREATE  (f%ice%work%fct_fluxes, f%ice%work%fct_plus, f%ice%work%fct_minus) 
    !$ACC ENTER DATA CREATE  (f%ice%work%eps11, f%ice%work%eps12, f%ice%work%eps22) 
    !$ACC ENTER DATA CREATE  (f%ice%work%sigma11, f%ice%work%sigma12, f%ice%work%sigma22) 
    !$ACC ENTER DATA CREATE  (f%ice%work%ice_strength, f%ice%stress_atmice_x, f%ice%stress_atmice_y) 
    !$ACC ENTER DATA COPYIN  (f%ice%thermo%rhosno, f%ice%thermo%rhoice, f%ice%thermo%inv_rhowat) 
    !$ACC ENTER DATA CREATE  (f%ice%srfoce_ssh, f%ice%pstar, f%ice%c_pressure) 
    !$ACC ENTER DATA CREATE  (f%ice%work%inv_areamass, f%ice%work%inv_mass, f%ice%uice_rhs, f%ice%vice_rhs) 
    !$ACC ENTER DATA CREATE  (f%ice%uice, f%ice%vice, f%ice%srfoce_u, f%ice%srfoce_v, f%ice%uice_old, f%ice%vice_old) 
    !$ACC ENTER DATA CREATE  (f%ice%data(1)%values, f%ice%data(2)%values, f%ice%data(3)%values) 
    !$ACC ENTER DATA CREATE  (f%ice%data(1)%valuesl, f%ice%data(2)%valuesl, f%ice%data(3)%valuesl) 
    !$ACC ENTER DATA CREATE  (f%ice%data(1)%dvalues, f%ice%data(2)%dvalues, f%ice%data(3)%dvalues) 
    !$ACC ENTER DATA CREATE  (f%ice%data(1)%values_rhs, f%ice%data(2)%values_rhs, f%ice%data(3)%values_rhs) 
    !$ACC ENTER DATA CREATE  (f%ice%data(1)%values_div_rhs, f%ice%data(2)%values_div_rhs, f%ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC ENTER DATA CREATE (f%ice%data(4)%values, f%ice%data(4)%valuesl, f%ice%data(4)%dvalues, f%ice%data(4)%values_rhs, f%ice%data(4)%values_div_rhs)
#endif
    !$ACC ENTER DATA COPYIN (f%dynamics)
    !$ACC ENTER DATA CREATE (f%dynamics%w, f%dynamics%w_e, f%dynamics%uv)
    !$ACC ENTER DATA CREATE (f%tracers%work%del_ttf)
    !$ACC ENTER DATA CREATE (f%tracers%data, f%tracers%work) 
    do tr_num=1, f%tracers%num_tracers
    !$ACC ENTER DATA CREATE (f%tracers%data(tr_num)%values, f%tracers%data(tr_num)%valuesAB)
    !$ACC ENTER DATA CREATE (f%tracers%data(tr_num)%tra_adv_ph, f%tracers%data(tr_num)%tra_adv_pv)
    end do
    !$ACC ENTER DATA CREATE (f%tracers%work%fct_ttf_min, f%tracers%work%fct_ttf_max, f%tracers%work%fct_plus, f%tracers%work%fct_minus) &
    !$ACC CREATE (f%tracers%work%adv_flux_hor, f%tracers%work%adv_flux_ver, f%tracers%work%fct_LO) &
    !$ACC CREATE (f%tracers%work%del_ttf_advvert, f%tracers%work%del_ttf_advhoriz, f%tracers%work%edge_up_dn_grad) &
    !$ACC CREATE (f%tracers%work%del_ttf)
  end subroutine


  subroutine fesom_runloop(current_nsteps)
    use fesom_main_storage_module
!   use openacc_lib
    integer, intent(in) :: current_nsteps 
    ! EO parameters
    integer n, nstart, ntotal, tr_num

    !=====================
    ! Time stepping
    !=====================

    ! --------------
    ! LA icebergs: 2023-05-17 
!YY: only when using icebergs?
    if (use_icebergs) then
    f%MPI_COMM_FESOM_IB = f%MPI_COMM_FESOM
    if (f%mype==0) then
!        write (*,*) 'ib_async_mode, initial omp_num_threads ', ib_async_mode, omp_get_num_threads()
        write (*,*) 'current_nsteps, steps_per_ib_step, icb_outfreq :', current_nsteps, steps_per_ib_step, icb_outfreq
    end if
    end if
    ! --------------
! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

    if (f%mype==0) write(*,*) 'FESOM start iteration before the barrier...'

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif 

    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)   
    if (f%mype==0) then

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

       write(*,*) 'FESOM start iteration after the barrier...'

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif 

       f%t0 = MPI_Wtime()
    endif

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

    if(f%mype==0) then
        write(*,*)
        print *, achar(27)//'[32m'  //'____________________________________________________________'//achar(27)//'[0m'
        print *, achar(27)//'[7;32m'//' --> FESOM STARTS TIME LOOP                                 '//achar(27)//'[0m'
    end if

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif

    !___MODEL TIME STEPPING LOOP________________________________________________
    if (use_global_tides) then
       call foreph_ini(yearnew, month, f%partit)
    end if
    
    nstart=f%from_nstep
    ntotal=f%from_nstep-1+current_nsteps

    do n=nstart, ntotal
        if (use_icebergs) then
                !n_ib         = n
                u_wind_ib    = u_wind
                v_wind_ib    = v_wind
                f%ice%uice_ib     = f%ice%uice
                f%ice%vice_ib     = f%ice%vice
        
        ! LA - this causes the blowup !
        !        f%ice%data(size(f%ice%data))      = f%ice%data(2)
        !        f%ice%data(size(f%ice%data)-1)    = f%ice%data(1)
        !!!!!!!!!!!!!!!!!!
        
        
        ! kh 08.03.21 support of different ocean ice and iceberg steps:
        ! if steps_per_ib_step is configured greater 1 then UV is modified via call oce_timestep_ale(n) -> call update_vel while
        ! the same asynchronous iceberg computation is still active
                f%dynamics%uv_ib     = f%dynamics%uv
        
        ! kh 15.03.21 support of different ocean ice and iceberg steps:
        ! if steps_per_ib_step is configured greater 1 then tr_arr is modified via call oce_timestep_ale(n) -> call solve_tracers_ale() while
        ! the same asynchronous iceberg computation is still active
                !tr_arr_ib    = tr_arr
                Tclim_ib     = f%tracers%data(1)%values
                Sclim_ib     = f%tracers%data(2)%values
        
        ! kh 15.03.21 support of different ocean ice and iceberg steps:
        ! if steps_per_ib_step is configured greater 1 then Tsurf and Ssurf might be changed while
        ! the same asynchronous iceberg computation is still active
                Tsurf_ib     = Tsurf
                Ssurf_ib     = Ssurf
        
        ! kh 18.03.21 support of different ocean ice and iceberg steps:
        ! if steps_per_ib_step is configured greater 1 then zbar_3d_n and eta_n might be changed while
        ! the same asynchronous iceberg computation is still active
                !zbar_3d_n_ib = zbar_3d_n
                f%mesh%Z_3d_n_ib     = f%mesh%Z_3d_n
                f%dynamics%eta_n_ib  = f%dynamics%eta_n
        
        ! kh 16.03.21 not modified during overlapping ocean/ice and iceberg computations
        !       coriolis_ib      = coriolis
        !       coriolis_node_ib = coriolis_node
        
        ! kh 02.02.21 check iceberg computations mode:
        ! ib_async_mode == 0: original sequential behavior for both ice sections (for testing purposes, creating reference results etc.)
        ! ib_async_mode == 1: OpenMP code active to overlapped computations in first (ocean ice) and second (icebergs) parallel section
        ! ib_async_mode == 2: OpenMP code active, but computations still serialized via spinlock (for testing purposes)
        
        ! -----------------------------------------------------------------------------------
        ! LA asyncronous coupling not included in this FESOM version, yet!!
        ! 
        end if        
        
        if (use_global_tides) then
           call foreph(f%partit, f%mesh)
        end if
        mstep = n

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif

        if (mod(n,logfile_outfreq)==0 .and. f%mype==0) then
            write(*,*) 'FESOM ======================================================='
!             write(*,*) 'FESOM step:',n,' day:', n*dt/24./3600.,
            write(*,*) 'FESOM step:',n,' day:', daynew,' year:',yearnew 
            write(*,*)
        end if

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif 

#if defined (__oifs) || defined (__oasis)
            seconds_til_now=INT(dt)*(n-1)
#endif
        call clock      
        !___compute horizontal velocity on nodes (originaly on elements)________
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_vel_nodes'//achar(27)//'[0m'
        call compute_vel_nodes(f%dynamics, f%partit, f%mesh)
        ! --------------
        ! LA icebergs: 2023-05-17 
        if (use_icebergs .and. mod(n - 1, steps_per_ib_step)==0) then
            if (f%mype==0) write(*,*) '*** step n=',n
            !t1_icb = MPI_Wtime()
            call iceberg_calculation(f%ice,f%mesh,f%partit,f%dynamics,n)
        end if
        ! --------------
        !___model sea-ice step__________________________________________________
        f%t1 = MPI_Wtime()
        if(use_ice) then
            !___compute fluxes from ocean to ice________________________________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ocean2ice(n)'//achar(27)//'[0m'
            call ocean2ice(f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
            
            !___compute update of atmospheric forcing____________________________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call update_atm_forcing(n)'//achar(27)//'[0m'
            f%t0_frc = MPI_Wtime()
            call update_atm_forcing(n, f%ice, f%tracers, f%dynamics, f%partit, f%mesh)
            f%t1_frc = MPI_Wtime()
            !___compute ice step________________________________________________
            if (f%ice%ice_steps_since_upd>=f%ice%ice_ave_steps-1) then
                f%ice%ice_update=.true.
                f%ice%ice_steps_since_upd = 0
            else
                f%ice%ice_update=.false.
                f%ice%ice_steps_since_upd=f%ice%ice_steps_since_upd+1
            endif
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call ice_timestep(n)'//achar(27)//'[0m'
            if (f%ice%ice_update) call ice_timestep(n, f%ice, f%partit, f%mesh)  

            
            ! LA commented for debugging
            ! --------------
            ! LA icebergs: 2023-05-17 
            if (use_icebergs .and. mod(n, steps_per_ib_step)==0.0) then
                call icb2fesom(f%mesh, f%partit, f%ice)
            end if
            ! --------------


            !___compute fluxes to the ocean: heat, freshwater, momentum_________
            if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_fluxes_mom...'//achar(27)//'[0m'
            call oce_fluxes_mom(f%ice, f%dynamics, f%partit, f%mesh) ! momentum only
            call oce_fluxes(f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        end if
        call before_oce_step(f%dynamics, f%tracers, f%partit, f%mesh) ! prepare the things if required
        f%t2 = MPI_Wtime()

        !___now recom____________________________________________________
#if defined (__recom)
! kh 29.03.22
#if defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

        if (f%mype==0 .and. n==1)  print *, achar(27)//'[46'  //'_____________________________________________________________'//achar(27)//'[0m'
        if (f%mype==0 .and. n==1)  print *, achar(27)//'[46;1m'//'     --> call REcoM                                         '//achar(27)//'[0m'

#if defined(__usetp)
        end if
#endif 

        f%t0_recom = MPI_Wtime()
        call recom(f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        f%t1_recom = MPI_Wtime()
#endif
        
        !___model ocean step____________________________________________________
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call oce_timestep_ale'//achar(27)//'[0m'
#if defined(__recom) && defined(__usetp)
        end if
#endif
        call oce_timestep_ale(n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)

        f%t3 = MPI_Wtime()
        !___compute energy diagnostics..._______________________________________
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call compute_diagnostics(1)'//achar(27)//'[0m'
#if defined(__recom) && defined(__usetp)
        end if
#endif
        call compute_diagnostics(1, f%dynamics, f%tracers, f%ice, f%partit, f%mesh)

        f%t4 = MPI_Wtime()
        !___prepare output______________________________________________________
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 
        if (flag_debug .and. f%mype==0)  print *, achar(27)//'[34m'//' --> call output (n)'//achar(27)//'[0m'
#if defined(__recom) && defined(__usetp)
        end if
#endif 
        call output (n, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)

        ! LA icebergs: 2023-05-17 
        if (use_icebergs .and. mod(n, steps_per_ib_step)==0.0) then
            call reset_ib_fluxes
        end if
        !--------------------------

        f%t5 = MPI_Wtime()
        call restart(n, nstart, f%total_nsteps, .false., f%which_readr, f%ice, f%dynamics, f%tracers, f%partit, f%mesh)
        f%t6 = MPI_Wtime()
        
        f%rtime_fullice       = f%rtime_fullice       + f%t2 - f%t1
        f%rtime_compute_diag  = f%rtime_compute_diag  + f%t4 - f%t3
        f%rtime_write_means   = f%rtime_write_means   + f%t5 - f%t4
        f%rtime_write_restart = f%rtime_write_restart + f%t6 - f%t5
        f%rtime_read_forcing  = f%rtime_read_forcing  + f%t1_frc - f%t0_frc
#if defined (__recom)
        f%rtime_compute_recom = f%rtime_compute_recom + f%t1_recom - f%t0_recom
#endif
    end do
!call cray_acc_set_debug_global_level(3)    
    f%from_nstep = f%from_nstep+current_nsteps
!call cray_acc_set_debug_global_level(0)    
!   write(0,*) 'f%from_nstep after the loop:', f%from_nstep    
  end subroutine


  subroutine fesom_finalize()
    use fesom_main_storage_module
#if defined(__MULTIO)
    use iom
    use mpp_io
#endif
    ! EO parameters
    real(kind=real32) :: mean_rtime(15), max_rtime(15), min_rtime(15)
    integer           :: tr_num
    integer :: i !YY
    ! --------------
    ! LA icebergs: 2023-05-17 
    if (use_icebergs) then
         call iceberg_out(f%partit)
    end if
    ! --------------
! kh 11.11.21 multi FESOM group loop parallelization
#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif 

    call finalize_output()
    call finalize_restart()

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif 

    !___FINISH MODEL RUN________________________________________________________

#if !defined (__usetp) 
! kh 11.11.21 multi FESOM group loop parallelization    
    call MPI_Barrier(f%MPI_COMM_FESOM, f%MPIERR)
#endif
#if defined(__recom) && defined (__usetp) 
! kh 11.11.21 list statistics for all fesom_groups 
! fesom groups are listed backwards, so info for the main fesom group 0 is at the end in the log
    do i = num_fesom_groups - 1, 0, -1

! kh 29.03.22 use a barrier to "sort" the output but the mpi output can still get a bit mixed up,
! because MPI does not define the handling of the order of the output lines
        call MPI_Barrier(f%MPI_COMM_FESOM_SAME_RANK_IN_GROUPS, f%MPIERR)

! kh 29.03.22 for the sake of output clarity produce output only for my_fesom_group == 0 for now
        if(i == f%my_fesom_group .and. f%my_fesom_group == 0) then
#endif

    !$ACC EXIT DATA DELETE (f%ice%delta_min, f%ice%Tevp_inv, f%ice%cd_oce_ice)
    !$ACC EXIT DATA DELETE (f%ice%work%fct_tmax, f%ice%work%fct_tmin)
    !$ACC EXIT DATA DELETE (f%ice%work%fct_fluxes, f%ice%work%fct_plus, f%ice%work%fct_minus)
    !$ACC EXIT DATA DELETE (f%ice%work%eps11, f%ice%work%eps12, f%ice%work%eps22)
    !$ACC EXIT DATA DELETE (f%ice%work%sigma11, f%ice%work%sigma12, f%ice%work%sigma22)
    !$ACC EXIT DATA DELETE (f%ice%work%ice_strength, f%ice%stress_atmice_x, f%ice%stress_atmice_y)
    !$ACC EXIT DATA DELETE (f%ice%thermo%rhosno, f%ice%thermo%rhoice, f%ice%thermo%inv_rhowat)
    !$ACC EXIT DATA DELETE (f%ice%srfoce_ssh, f%ice%pstar, f%ice%c_pressure)
    !$ACC EXIT DATA DELETE (f%ice%work%inv_areamass, f%ice%work%inv_mass, f%ice%uice_rhs, f%ice%vice_rhs)
    !$ACC EXIT DATA DELETE (f%ice%uice, f%ice%vice, f%ice%srfoce_u, f%ice%srfoce_v, f%ice%uice_old, f%ice%vice_old)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%values, f%ice%data(2)%values, f%ice%data(3)%values)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%valuesl, f%ice%data(2)%valuesl, f%ice%data(3)%valuesl)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%dvalues, f%ice%data(2)%dvalues, f%ice%data(3)%dvalues)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%values_rhs, f%ice%data(2)%values_rhs, f%ice%data(3)%values_rhs)
    !$ACC EXIT DATA DELETE (f%ice%data(1)%values_div_rhs, f%ice%data(2)%values_div_rhs, f%ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC EXIT DATA DELETE (f%ice%data(4)%values, f%ice%data(4)%valuesl, f%ice%data(4)%dvalues, f%ice%data(4)%values_rhs, f%ice%data(4)%values_div_rhs)
#endif
    !$ACC EXIT DATA DELETE (f%ice%data, f%ice%work, f%ice%work%fct_massmatrix)
    !$ACC EXIT DATA DELETE (f%ice)
    do tr_num=1, f%tracers%num_tracers
    !$ACC EXIT DATA DELETE (f%tracers%data(tr_num)%values, f%tracers%data(tr_num)%valuesAB)
    end do
    !$ACC EXIT DATA DELETE (f%tracers%work%fct_ttf_min, f%tracers%work%fct_ttf_max, f%tracers%work%fct_plus, f%tracers%work%fct_minus)
    !$ACC EXIT DATA DELETE (f%tracers%work%adv_flux_hor, f%tracers%work%adv_flux_ver, f%tracers%work%fct_LO)
    !$ACC EXIT DATA DELETE (f%tracers%work%del_ttf_advvert, f%tracers%work%del_ttf_advhoriz, f%tracers%work%edge_up_dn_grad)
    !$ACC EXIT DATA DELETE (f%tracers%work%del_ttf)
    !$ACC EXIT DATA DELETE (f%tracers%data, f%tracers%work)
    !$ACC EXIT DATA DELETE (f%dynamics%w, f%dynamics%w_e, f%dynamics%uv)
    !$ACC EXIT DATA DELETE (f%dynamics, f%tracers)

    !delete mesh and partit data.
    !$ACC EXIT DATA DELETE (f%mesh%coriolis_node, f%mesh%nn_num, f%mesh%nn_pos) 
    !$ACC EXIT DATA DELETE (f%mesh%ssh_stiff, f%mesh%ssh_stiff%rowptr) 
    !$ACC EXIT DATA DELETE (f%mesh%gradient_sca, f%mesh%metric_factor, f%mesh%elem_area, f%mesh%area, f%mesh%edge2D_in) 
    !$ACC EXIT DATA DELETE (f%mesh%elem2D_nodes, f%mesh%ulevels, f%mesh%ulevels_nod2d, f%mesh%edges, f%mesh%edge_tri) 
    !$ACC EXIT DATA DELETE (f%mesh%helem, f%mesh%elem_cos, f%mesh%edge_cross_dxdy, f%mesh%elem2d_nodes, f%mesh%nl) 
    !$ACC EXIT DATA DELETE (f%mesh%nlevels_nod2D, f%mesh%nod_in_elem2D, f%mesh%nod_in_elem2D_num) 
    !$ACC EXIT DATA DELETE (f%mesh%edge_dxdy, f%mesh%nlevels, f%mesh%hnode, f%mesh%hnode_new, f%mesh%ulevels_nod2D_max) 
    !$ACC EXIT DATA DELETE (f%mesh%zbar_3d_n, f%mesh%z_3d_n, f%mesh%areasvol, f%mesh%nlevels_nod2D_min) 
    !$ACC EXIT DATA DELETE (f%partit%eDim_nod2D, f%partit%myDim_edge2D) 
    !$ACC EXIT DATA DELETE (f%partit%myDim_elem2D, f%partit%myDim_nod2D, f%partit%myList_edge2D) 
    !$ACC EXIT DATA DELETE (f%mesh, f%partit, f)
    if (f%mype==0) then
       f%t1 = MPI_Wtime()
       f%runtime_alltimesteps = real(f%t1-f%t0,real32)
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
    mean_rtime(10) = f%rtime_fullice - f%rtime_read_forcing 
    mean_rtime(11) = f%rtime_compute_diag
    mean_rtime(12) = f%rtime_write_means
    mean_rtime(13) = f%rtime_write_restart
    mean_rtime(14) = f%rtime_read_forcing
#if defined (__recom)
    mean_rtime(15) = f%rtime_compute_recom
#endif
    max_rtime(1:14) = mean_rtime(1:14)
    min_rtime(1:14) = mean_rtime(1:14)
#if defined (__recom)
    max_rtime(15) = mean_rtime(15)
    min_rtime(15) = mean_rtime(15)
    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime(15), 1, MPI_REAL, MPI_SUM, f%MPI_COMM_FESOM, f%MPIerr)
    mean_rtime(15) = mean_rtime(15) / real(f%npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime(15),  1, MPI_REAL, MPI_MAX, f%MPI_COMM_FESOM, f%MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime(15),  1, MPI_REAL, MPI_MIN, f%MPI_COMM_FESOM, f%MPIerr)
#endif

    call MPI_AllREDUCE(MPI_IN_PLACE, mean_rtime, 14, MPI_REAL, MPI_SUM, f%MPI_COMM_FESOM, f%MPIerr)
    mean_rtime(1:14) = mean_rtime(1:14) / real(f%npes,real32)
    call MPI_AllREDUCE(MPI_IN_PLACE, max_rtime,  14, MPI_REAL, MPI_MAX, f%MPI_COMM_FESOM, f%MPIerr)
    call MPI_AllREDUCE(MPI_IN_PLACE, min_rtime,  14, MPI_REAL, MPI_MIN, f%MPI_COMM_FESOM, f%MPIerr)
    
#if defined (__oifs) 
    ! OpenIFS coupled version has to call oasis_terminate through par_ex
    call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype)
#endif

#if defined(__recom) && defined (__usetp)
        end if
    end do ! i = num_fesom_groups - 1, 0, -1
#endif

#if defined(__MULTIO) && !defined(__ifsinterface) && !defined(__oasis)
   call mpp_stop
#endif
    if(f%fesom_did_mpi_init) call par_ex(f%partit%MPI_COMM_FESOM, f%partit%mype) ! finalize MPI before FESOM prints its stats block, otherwise there is sometimes output from other processes from an earlier time in the programm AFTER the starts block (with parastationMPI)

#if defined(__recom) && defined(__usetp)
        if (f%my_fesom_group==0) then
#endif

    if (f%mype==0) then
        41 format (a35,a10,2a15) !Format for table heading
        42 format (a30,3f15.4)   !Format for table content

        print 41, '___MODEL RUNTIME per task [seconds]','_____mean_','___________min_', '___________max_'
        print 42, '  runtime ocean:              ',    mean_rtime(1),     min_rtime(1),      max_rtime(1)
        print 42, '    > runtime oce. mix,pres. :',    mean_rtime(2),     min_rtime(2),      max_rtime(2)
        print 42, '    > runtime oce. dyn. u,v,w:',    mean_rtime(3),     min_rtime(3),      max_rtime(3)
        print 42, '    > runtime oce. dyn. ssh  :',    mean_rtime(4),     min_rtime(4),      max_rtime(4)
        print 42, '    > runtime oce. solve ssh :',    mean_rtime(5),     min_rtime(5),      max_rtime(5)
        print 42, '    > runtime oce. GM/Redi   :',    mean_rtime(6),     min_rtime(6),      max_rtime(6)
        print 42, '    > runtime oce. tracer    :',    mean_rtime(7),     min_rtime(7),      max_rtime(7)
        print 42, '  runtime ice  :              ',    mean_rtime(10),    min_rtime(10),     max_rtime(10)
        print 42, '    > runtime ice step :      ',    mean_rtime(8),     min_rtime(8),      max_rtime(8)
        print 42, '  runtime diag:               ',    mean_rtime(11),    min_rtime(11),     max_rtime(11)
        print 42, '  runtime output:             ',    mean_rtime(12),    min_rtime(12),     max_rtime(12)
        print 42, '  runtime restart:            ',    mean_rtime(13),    min_rtime(13),     max_rtime(13)
        print 42, '  runtime forcing:            ',    mean_rtime(14),    min_rtime(14),     max_rtime(14)
        print 42, '  runtime total (ice+oce):    ',    mean_rtime(9),     min_rtime(9),      max_rtime(9)
#if defined (__recom)
        print 42, '  runtime recom:              ',    mean_rtime(15),    min_rtime(15),     max_rtime(15)
#endif

        43 format (a33,i15)        !Format Ncores
        44 format (a33,i15)        !Format OMP threads
        45 format (a33,f15.4,a4)   !Format runtime

        write(*,*)
        write(*,*) '======================================================'
        write(*,*) '================ BENCHMARK RUNTIME ==================='
        print 43, '    Number of cores :            ',f%npes
#if defined(_OPENMP)
        print 44, '    Max OpenMP threads :         ',OMP_GET_MAX_THREADS()
#endif
        print 45, '    Runtime for all timesteps :  ',f%runtime_alltimesteps,' sec'
        write(*,*) '======================================================'
        write(*,*)
    end if    

! kh 29.03.22
#if defined(__recom) && defined(__usetp)
        end if
#endif 

!   call clock_finish  
  end subroutine

end module
