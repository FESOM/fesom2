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
        integer         intent(in)            :: istep
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
    use ice_array_setup_interface
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
    integer,        intent(in)            :: step
    type(t_ice),    intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh),   intent(in),    target :: mesh
    !___________________________________________________________________________
    integer                               :: i
    REAL(kind=WP)                         :: t0,t1, t2, t3
#if defined (__icepack)
    real(kind=WP)                         :: time_evp, time_advec, time_therm
#endif
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
#if defined (__oifs) || defined (__ifsinterface)
    real(kind=WP), dimension(:), pointer  :: ice_temp, a_ice
#endif     
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    u_ice    => ice%uice(:)
    v_ice    => ice%vice(:)
#if defined (__oifs) || defined (__ifsinterface)
    a_ice    => ice%data(1)%values(:)
    ice_temp => ice%data(4)%values(:)
#endif     
    !___________________________________________________________________________
    t0=MPI_Wtime()
#if defined (__icepack)
    call step_icepack(ice, mesh, time_evp, time_advec, time_therm) ! EVP, advection and thermodynamic parts    
#else     
    
    !___________________________________________________________________________
    ! ===== Dynamics
    SELECT CASE (ice%whichEVP)
    CASE (0)
        if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'  
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
    do i=1,myDim_nod2D+eDim_nod2D
        ice_temp(i) = ice_temp(i)*a_ice(i)
    end do
#endif /* (__oifs) */
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_TG_rhs_div...'//achar(27)//'[0m'
    call ice_TG_rhs_div    (ice, partit, mesh)  
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m' 
    call ice_fct_solve     (ice, partit, mesh)
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
    call ice_update_for_div(ice, partit, mesh)
    
#if defined (__oifs) || defined (__ifsinterface)
    do i=1,myDim_nod2D+eDim_nod2D
        if (a_ice(i)>0.0_WP) ice_temp(i) = ice_temp(i)/a_ice(i)
    end do
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
    do i=1,myDim_nod2D+eDim_nod2D
        if ( ( U_ice(i)/=0.0_WP .and. mesh%ulevels_nod2d(i)>1) .or. (V_ice(i)/=0.0_WP .and. mesh%ulevels_nod2d(i)>1) ) then
            write(*,*) " --> found cavity velocity /= 0.0_WP , ", mype
            write(*,*) " ulevels_nod2d(n) = ", mesh%ulevels_nod2d(i)
            write(*,*) " U_ice(n) = ", U_ice(i)
            write(*,*) " V_ice(n) = ", V_ice(i)
            write(*,*)
        end if 
    end do
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
    implicit none
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(in)   , target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: i
    character(MAX_PATH)                   :: filename
    real(kind=WP), external               :: TFrez  ! Sea water freeze temperature.
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:), pointer  :: a_ice, m_ice, m_snow
    real(kind=WP), dimension(:), pointer  :: u_ice, v_ice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
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
    m_snow=0._WP
    if(mype==0) write(*,*) 'initialize the sea ice'
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
end subroutine ice_initial_state
