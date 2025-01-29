module ice_initial_state_interface
    interface
        subroutine ice_initial_state(ice, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_TRACER
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_tracer), intent(in)   , target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module ice_setup_interface
    interface
        subroutine ice_setup(ice, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_TRACER
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        type(t_ice)   , intent(inout), target :: ice
        type(t_tracer), intent(in)   , target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module ice_timestep_interface
    interface
        subroutine ice_timestep(istep, ice, partit, mesh)
        USE MOD_ICE
        USE MOD_PARTIT
        USE MOD_PARSUP
        USE MOD_MESH
        integer       , intent(in)            :: istep
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

!
!_______________________________________________________________________________
! ice initialization + array allocation + time stepping
subroutine ice_setup(ice, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_CONFIG
    use ice_initial_state_interface
    use ice_fct_interfaces
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(in)   , target :: tracers
    type(t_mesh)  , intent(in)   , target :: mesh
    type(t_partit), intent(inout), target :: partit

    !___________________________________________________________________________
    ! initialise ice derived type
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_init'//achar(27)//'[0m'
    call ice_init(ice, partit, mesh)

    !___________________________________________________________________________
    ! DO not change
    ice%ice_dt   = real(ice%ice_ave_steps,WP)*dt
    ! ice_dt=dt
    ice%Tevp_inv = 3.0_WP/ice%ice_dt
    ! This is combination it always enters
    ice%Clim_evp = ice%Clim_evp*(ice%evp_rheol_steps/ice%ice_dt)**2/ice%Tevp_inv

    !___________________________________________________________________________
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_init'//achar(27)//'[0m'
    call ice_mass_matrix_fill(ice, partit, mesh)

    !___________________________________________________________________________
    ! Initialization routine, user input is required
    !call ice_init_fields_test
    if (flag_debug .and. partit%mype==0)  print *, achar(27)//'[36m'//'     --> call ice_initial_state'//achar(27)//'[0m'
    call ice_initial_state(ice, tracers, partit, mesh)   ! Use it unless running test example

    if(partit%mype==0) write(*,*) 'Ice is initialized'
end subroutine ice_setup
!
!
!_______________________________________________________________________________
! Sea ice model step
subroutine ice_timestep(step, ice, partit, mesh)
    USE MOD_ICE
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_param
    use g_CONFIG
    use ice_EVPdynamics_interface
    use ice_maEVPdynamics_interface
    use ice_fct_interfaces
    use ice_thermodynamics_interfaces
    use cavity_interfaces
#if defined (__icepack)
    use icedrv_main,   only: step_icepack
#endif
    implicit none
    integer       , intent(in)            :: step
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: i
    REAL(kind=WP)                         :: t0,t1, t2, t3
#if defined (__icepack)
    real(kind=WP)                         :: time_evp, time_advec, time_therm
#endif
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    !LA 2023-03-08
    real(kind=WP), dimension(:), pointer  :: u_ice_ib, v_ice_ib
#if defined (__oifs) || defined (__ifsinterface)
    real(kind=WP), dimension(:), pointer  :: a_ice, ice_temp
    !LA 2023-03-08
    real(kind=WP), dimension(:), pointer  :: a_ice_ib
#endif
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

!---------------------------------------------
! LA: 2023-01-31 add asynchronous icebergs
!    u_ice    => ice%uice(:)
!    v_ice    => ice%vice(:)
! kh 19.02.21
  if (ib_async_mode == 0) then
      u_ice    => ice%uice(:)
      v_ice    => ice%vice(:)
      !allocate(u_ice(n_size), v_ice(n_size))
      !allocate(u_ice_ib(n_size), v_ice_ib(n_size))
  else
!$omp parallel sections num_threads(2)
! kh 19.02.21 support "first touch" idea
!$omp section
      u_ice    => ice%uice(:)
      v_ice    => ice%vice(:)
      !allocate(u_ice(n_size), v_ice(n_size))
      u_ice    = 0._WP
      v_ice    = 0._WP
      !do i = 1, n_size
      !    u_ice(i) = 0._WP
      !    v_ice(i) = 0._WP
      !end do
!$omp section
      u_ice_ib => ice%uice_ib(:)
      v_ice_ib => ice%vice_ib(:)
      !allocate(u_ice_ib(n_size), v_ice_ib(n_size))
      u_ice_ib = 0._WP
      v_ice_ib = 0._WP
      !do i = 1, n_size
      !    u_ice_ib(i) = 0._WP
      !    v_ice_ib(i) = 0._WP
      !end do
!$omp end parallel sections
  end if
!---------------------------------------------
#if defined (__oifs) || defined (__ifsinterface)
    a_ice    => ice%data(1)%values(:)    
    ice_temp => ice%data(4)%values(:)
#endif
    !___________________________________________________________________________
    t0=MPI_Wtime()
#if defined (__icepack)
    call step_icepack(ice, mesh, time_evp, time_advec, time_therm) ! EVP, advection and thermodynamic parts
#else

    !$ACC UPDATE DEVICE (ice%work%fct_massmatrix) &
    !$ACC DEVICE (ice%delta_min, ice%Tevp_inv, ice%cd_oce_ice) &
    !$ACC DEVICE (ice%work%fct_tmax, ice%work%fct_tmin) &
    !$ACC DEVICE (ice%work%fct_fluxes, ice%work%fct_plus, ice%work%fct_minus) &
    !$ACC DEVICE (ice%work%eps11, ice%work%eps12, ice%work%eps22) &
    !$ACC DEVICE (ice%work%sigma11, ice%work%sigma12, ice%work%sigma22) &
    !$ACC DEVICE (ice%work%ice_strength, ice%stress_atmice_x, ice%stress_atmice_y) &
    !$ACC DEVICE (ice%thermo%rhosno, ice%thermo%rhoice, ice%thermo%inv_rhowat) &
    !$ACC DEVICE (ice%srfoce_ssh, ice%pstar, ice%c_pressure) &
    !$ACC DEVICE (ice%work%inv_areamass, ice%work%inv_mass, ice%uice_rhs, ice%vice_rhs) &
    !$ACC DEVICE (ice%uice, ice%vice, ice%srfoce_u, ice%srfoce_v, ice%uice_old, ice%vice_old) &
    !$ACC DEVICE (ice%data(1)%values, ice%data(2)%values, ice%data(3)%values) &
    !$ACC DEVICE (ice%data(1)%valuesl, ice%data(2)%valuesl, ice%data(3)%valuesl) &
    !$ACC DEVICE (ice%data(1)%dvalues, ice%data(2)%dvalues, ice%data(3)%dvalues) &
    !$ACC DEVICE (ice%data(1)%values_rhs, ice%data(2)%values_rhs, ice%data(3)%values_rhs) &
    !$ACC DEVICE (ice%data(1)%values_div_rhs, ice%data(2)%values_div_rhs, ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC UPDATE DEVICE (ice%data(4)%values, ice%data(4)%valuesl, ice%data(4)%dvalues, ice%data(4)%values_rhs, ice%data(4)%values_div_rhs)
#endif
    !___________________________________________________________________________
    ! ===== Dynamics

    SELECT CASE (ice%whichEVP)
    CASE (0)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'
#if defined(_CRAYFTN)
	!dir$ noinline
#endif
        call EVPdynamics  (ice, partit, mesh)
    CASE (1)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_m...'//achar(27)//'[0m'
        call EVPdynamics_m(ice, partit, mesh)
    CASE (2)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics_a...'//achar(27)//'[0m'
        call EVPdynamics_a(ice, partit, mesh)
    CASE DEFAULT
        if (mype==0) write(*,*) 'a non existing EVP scheme specified!'
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    END SELECT

    if (use_cavity) call cavity_ice_clean_vel(ice, partit, mesh)
    t1=MPI_Wtime()

    !___________________________________________________________________________
    ! ===== Advection part
    ! old FCT routines
    ! call ice_TG_rhs
    ! call ice_fct_solve
    ! call cut_off
    ! new FCT routines from Sergey Danilov 08.05.2018
#if defined (__oifs) || defined (__ifsinterface)
#ifndef ENABLE_OPENACC
!$OMP PARALLEL DO
#else
!$ACC parallel loop present(ice_temp, a_ice)
#endif
    do i=1,myDim_nod2D+eDim_nod2D
        ice_temp(i) = ice_temp(i)*a_ice(i)
    end do
#ifndef ENABLE_OPENACC
!$OMP END PARALLEL DO
#else
!$ACC END parallel loop
#endif
#endif /* (__oifs) */
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_TG_rhs_div...'//achar(27)//'[0m'
    call ice_TG_rhs    (ice, partit, mesh)

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m'
    call ice_fct_solve     (ice, partit, mesh)

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
!   call ice_update_for_div(ice, partit, mesh)

    !$ACC UPDATE HOST (ice%work%fct_massmatrix) &
    !$ACC HOST (ice%delta_min, ice%Tevp_inv, ice%cd_oce_ice) &
    !$ACC HOST (ice%work%fct_tmax, ice%work%fct_tmin) &
    !$ACC HOST (ice%work%fct_fluxes, ice%work%fct_plus, ice%work%fct_minus) &
    !$ACC HOST (ice%work%eps11, ice%work%eps12, ice%work%eps22) &
    !$ACC HOST (ice%work%sigma11, ice%work%sigma12, ice%work%sigma22) &
    !$ACC HOST (ice%work%ice_strength, ice%stress_atmice_x, ice%stress_atmice_y) &
    !$ACC HOST (ice%thermo%rhosno, ice%thermo%rhoice, ice%thermo%inv_rhowat) &
    !$ACC HOST (ice%srfoce_ssh, ice%pstar, ice%c_pressure) &
    !$ACC HOST (ice%work%inv_areamass, ice%work%inv_mass, ice%uice_rhs, ice%vice_rhs) &
    !$ACC HOST (ice%uice, ice%vice, ice%srfoce_u, ice%srfoce_v, ice%uice_old, ice%vice_old) &
    !$ACC HOST (ice%data(1)%values, ice%data(2)%values, ice%data(3)%values) &
    !$ACC HOST (ice%data(1)%valuesl, ice%data(2)%valuesl, ice%data(3)%valuesl) &
    !$ACC HOST (ice%data(1)%dvalues, ice%data(2)%dvalues, ice%data(3)%dvalues) &
    !$ACC HOST (ice%data(1)%values_rhs, ice%data(2)%values_rhs, ice%data(3)%values_rhs) &
    !$ACC HOST (ice%data(1)%values_div_rhs, ice%data(2)%values_div_rhs, ice%data(3)%values_div_rhs)
#if defined (__oifs) || defined (__ifsinterface)
    !$ACC UPDATE HOST (ice%data(4)%values, ice%data(4)%valuesl, ice%data(4)%dvalues, ice%data(4)%values_rhs, ice%data(4)%values_div_rhs)
#endif

#if defined (__oifs) || defined (__ifsinterface)
!$OMP PARALLEL DO
    do i=1,myDim_nod2D+eDim_nod2D
        if (a_ice(i)>0.0_WP) ice_temp(i) = ice_temp(i)/max(a_ice(i), 1.e-6_WP)
    end do
!$OMP END PARALLEL DO
#endif /* (__oifs) */

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call cut_off...'//achar(27)//'[0m'
    call cut_off(ice, partit, mesh)

    if (use_cavity) call cavity_ice_clean_ma(ice, partit, mesh)
    t2=MPI_Wtime()

    !___________________________________________________________________________
    ! ===== Thermodynamic part
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call thermodynamics...'//achar(27)//'[0m'
    call thermodynamics(ice, partit, mesh)
#endif /* (__icepack) */


    !___________________________________________________________________________
!$OMP PARALLEL DO
    do i=1,myDim_nod2D+eDim_nod2D
        ice%h_ice(i) =ice%data(2)%values(i)/max(ice%data(1)%values(i), 1.e-3)
        ice%h_snow(i)=ice%data(3)%values(i)/max(ice%data(1)%values(i), 1.e-3)
        if ( ( U_ice(i)/=0.0_WP .and. mesh%ulevels_nod2d(i)>1) .or. (V_ice(i)/=0.0_WP .and. mesh%ulevels_nod2d(i)>1) ) then
            write(*,*) " --> found cavity velocity /= 0.0_WP , ", mype
            write(*,*) " ulevels_nod2d(n) = ", mesh%ulevels_nod2d(i)
            write(*,*) " U_ice(n) = ", U_ice(i)
            write(*,*) " V_ice(n) = ", V_ice(i)
            write(*,*)
        end if
    end do
!$OMP END PARALLEL DO
    t3=MPI_Wtime()
    rtime_ice = rtime_ice + (t3-t0)
    rtime_tot = rtime_tot + (t3-t0)
    if(mod(step,logfile_outfreq)==0 .and. mype==0) then
        write(*,*) '___ICE STEP EXECUTION TIMES____________________________'
#if defined (__icepack)
        write(*,"(A, ES10.3)") '	Ice Dyn.        :', time_evp
                write(*,"(A, ES10.3)") '        Ice Advect.     :', time_advec
                write(*,"(A, ES10.3)") '        Ice Thermodyn.  :', time_therm
#else
        write(*,"(A, ES10.3)") '	Ice Dyn.        :', t1-t0
        write(*,"(A, ES10.3)") '	Ice Advect.     :', t2-t1
        write(*,"(A, ES10.3)") '	Ice Thermodyn.  :', t3-t2
#endif /* (__icepack) */
        write(*,*) '   _______________________________'
        write(*,"(A, ES10.3)") '	Ice TOTAL       :', t3-t0
        write(*,*)
     endif
end subroutine ice_timestep
!
!
!_______________________________________________________________________________
! sets inital values or reads restart file for ice model
subroutine ice_initial_state(ice, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    USE MOD_MESH
    use o_PARAM
    use o_arrays
    use g_CONFIG
    USE g_read_other_NetCDF, only: read_other_NetCDF
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(in)   , target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: i
    character(MAX_PATH)                   :: filename
    real(kind=WP), external               :: TFrez  ! Sea water freeze temperature.
!============== namelistatmdata variables ================
   integer, save                                :: nm_ic_unit     = 107 ! unit to open namelist file
   integer                                      :: iost                 !I/O status
   integer, parameter                           :: ic_max=10
   logical                                      :: ic_cyclic=.true.
   integer,             save                    :: n_ic2d
   integer,             save, dimension(ic_max) :: idlist
   character(MAX_PATH), save, dimension(ic_max) :: filelist
   logical                                      :: ini_ice_from_file=.false.
   character(50),       save, dimension(ic_max) :: varlist
   integer                                      :: current_tracer
   namelist / tracer_init2d / n_ic2d, idlist, filelist, varlist, ini_ice_from_file

    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
    !LA 2023-03-07
    real(kind=WP), dimension(:), pointer  :: a_ice_ib, m_ice_ib
    real(kind=WP), dimension(:), pointer  :: u_ice_ib, v_ice_ib
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"

! LA: 2023-01-31 add asynchronous icebergs
!---------------------------------------------
m_snow       => ice%data(3)%values(:)    
m_snow=0._WP
if (.not.use_icebergs) then
    u_ice        => ice%uice(:)
    v_ice        => ice%vice(:)
    a_ice        => ice%data(1)%values(:)
    m_ice        => ice%data(2)%values(:)
    m_snow       => ice%data(3)%values(:)
    !___________________________________________________________________________
    m_ice =0._WP
    a_ice =0._WP
    u_ice =0._WP
    v_ice =0._WP
else
  if (ib_async_mode == 0) then
    u_ice        => ice%uice(:)
    v_ice        => ice%vice(:)
    u_ice_ib     => ice%uice_ib(:)
    v_ice_ib     => ice%vice_ib(:)
    a_ice        => ice%data(1)%values(:)
    m_ice        => ice%data(2)%values(:)
    a_ice_ib     => ice%data(size(ice%data)-1)%values(:)
    m_ice_ib     => ice%data(size(ice%data))%values(:)
    !allocate(m_ice(n_size), a_ice(n_size))
    !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
    m_ice        = 0._WP
    a_ice        = 0._WP
    u_ice        = 0._WP
    v_ice        = 0._WP
    u_ice_ib     = 0._WP
    v_ice_ib     = 0._WP
    m_ice_ib     = 0._WP
    a_ice_ib     = 0._WP
  else
! kh 19.02.21 support "first touch" idea
!$omp parallel sections num_threads(2)
!$omp section
      !allocate(m_ice(n_size), a_ice(n_size))
      !do i = 1, n_size
      !    m_ice(i) = 0._WP
      !    a_ice(i) = 0._WP
      !end do
      u_ice        => ice%uice(:)
      v_ice        => ice%vice(:)
      a_ice        => ice%data(1)%values(:)
      m_ice        => ice%data(2)%values(:)
      !allocate(m_ice(n_size), a_ice(n_size))
      !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
      m_ice        = 0._WP
      a_ice        = 0._WP
      u_ice        = 0._WP
      v_ice        = 0._WP
!$omp section
      !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
      !do i = 1, n_size
      !    m_ice_ib(i) = 0._WP
      !    a_ice_ib(i) = 0._WP
      !end do
      u_ice_ib     => ice%uice_ib(:)
      v_ice_ib     => ice%vice_ib(:)
      a_ice_ib     => ice%data(size(ice%data)-1)%values(:)
      m_ice_ib     => ice%data(size(ice%data))%values(:)
      !allocate(m_ice(n_size), a_ice(n_size))
      !allocate(m_ice_ib(n_size), a_ice_ib(n_size))
      u_ice_ib     = 0._WP
      v_ice_ib     = 0._WP
      m_ice_ib     = 0._WP
      a_ice_ib     = 0._WP
!$omp end parallel sections
  end if
end if
! LA: 2023-01-31 add asynchronous icebergs
!---------------------------------------------



    !___________________________________________________________________________
    ! OPEN and read namelist for I/O
    open( unit=nm_ic_unit, file='namelist.tra', form='formatted', access='sequential', status='old', iostat=iost )
    if (iost == 0) then
        if (mype==0) WRITE(*,*) '     file   : ', 'namelist.tra',' open ok'
        else
        if (mype==0) WRITE(*,*) 'ERROR: --> bad opening file   : ', 'namelist.tra',' ; iostat=',iost
        call par_ex(partit%MPI_COMM_FESOM, partit%mype)
        stop
    end if
    read(nm_ic_unit, nml=tracer_init2d,   iostat=iost)
    close(nm_ic_unit)

    if (.not. ini_ice_from_file) then
    if(mype==0) write(*,*) 'initialize the sea ice: cold start'
    !___________________________________________________________________________
    do i=1,myDim_nod2D+eDim_nod2D
        !_______________________________________________________________________
        if (ulevels_nod2d(i)>1) then
            !!PS m_ice(i)  = 1.0e15_WP
            !!PS m_snow(i) = 0.1e15_WP
            cycle ! --> if cavity, no sea ice, no initial state
        endif

        !_______________________________________________________________________
        if (tracers%data(1)%values(1,i)< 0.0_WP) then
            if (geo_coord_nod2D(2,i)>0._WP) then
                m_ice(i) = 1.0_WP
                m_snow(i)= 0.1_WP
            else
                m_ice(i) = 2.0_WP
                m_snow(i)= 0.5_WP
            end if
            a_ice(i) = 0.9_WP
            u_ice(i) = 0.0_WP
            v_ice(i) = 0.0_WP
        endif
    enddo
    else
    if (mype==0) write(*,*) 'initialize the sea ice: from file'
       DO i=1, n_ic2d
          DO current_tracer=1, ice%num_itracers
             IF (ice%data(current_tracer)%ID==idlist(i)) then
                IF (mype==0) then
                   write(*,*) 'reading 2D variable  : ', trim(varlist(i)), ' into 2D tracer ID=', current_tracer
                   write(*,*) 'from ',                   trim(filelist(i))
                END IF
             call read_other_NetCDF(trim(ClimateDataPath)//trim(filelist(i)), varlist(i),  1, ice%data(current_tracer)%values(:), .false., .true., partit, mesh)
             EXIT
             END IF
          END DO
       END DO
    end if
end subroutine ice_initial_state
