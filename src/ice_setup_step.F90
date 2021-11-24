module ice_setup_interface
    interface
        subroutine ice_setup(ice, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_TRACER
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_ice)   , intent(inout), target :: ice
        type(t_tracer), intent(in)   , target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module ice_initial_state_interface
    interface
        subroutine ice_initial_state(ice, tracers, partit, mesh)
        USE MOD_ICE
        USE MOD_TRACER
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        type(t_ice)   , intent(inout), target :: ice
        type(t_tracer), intent(in)   , target :: tracers
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module

module ice_timestep_interface
    interface
        subroutine ice_timestep(step, ice, partit, mesh)
        USE MOD_ICE
        USE MOD_MESH
        USE MOD_PARTIT
        USE MOD_PARSUP
        integer                               :: step
        type(t_ice)   , intent(inout), target :: ice
        type(t_partit), intent(inout), target :: partit
        type(t_mesh)  , intent(in)   , target :: mesh
        end subroutine
    end interface
end module
!
!
!_______________________________________________________________________________
! ice initialization + array allocation + time stepping
subroutine ice_setup(ice, tracers, partit, mesh)
    USE MOD_ICE
    USE MOD_TRACER
    USE MOD_PARTIT
    USE MOD_PARSUP
    use MOD_MESH
    use o_param
    use i_param
    use i_arrays
    use g_CONFIG
    use ice_init_interface
    use ice_fct_interfaces
    use ice_initial_state_interface
    implicit none 
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(in)   , target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    ! pointer on necessary derived types
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    !___________________________________________________________________________
    ! initialise ice derived type 
    call ice_init(ice, partit, mesh)
    
    !___________________________________________________________________________
    ! DO not change: ice_dt=dt
    ice%ice_dt=real(ice%ice_ave_steps,WP)*dt
    ice%Tevp_inv=3.0_WP/ice%ice_dt 
    ! This is combination it always enters
    ice%Clim_evp=ice%Clim_evp*(ice%evp_rheol_steps/ice%ice_dt)**2/ice%Tevp_inv  
    
    !___________________________________________________________________________
    ! Fill in  the mass matrix  
    call ice_fct_init(ice, partit, mesh)
!     call ice_mass_matrix_fill(ice, partit, mesh)
    
    !___________________________________________________________________________
    ! Initialization routine, user input is required 
    !call ice_init_fields_test
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
    use i_arrays
    use o_param, only: WP
    use g_CONFIG
!     use i_PARAM, only: whichEVP
    use ice_EVP_interfaces
    use ice_maEVP_interfaces
    use ice_TG_rhs_div_interfaces
    use ice_update_for_div_interface
    use ice_fct_solve_interface
    use ice_thermodynamics_interfaces
    use cavity_ice_clean_interface
#if defined (__icepack)
    use icedrv_main,   only: step_icepack 
#endif
    implicit none 
    integer                               :: step
    type(t_ice)   , intent(inout), target :: ice
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: node
    real(kind=WP)                         :: t0, t1, t2, t3    
#if defined (__icepack)
    real(kind=WP)                         :: time_evp, time_advec, time_therm
#endif
    !___________________________________________________________________________
    ! pointer on necessary derived types
#if defined (__oifs)    
    real(kind=WP), dimension(:), pointer :: ice_temp, a_ice
#endif /* (__oifs) */    
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
#if defined (__oifs) 
    a_ice   => ice%data(1)%values(:)
    ice_temp=> ice%data(4)%values(:)
#endif /* (__oifs) */    

    !___________________________________________________________________________
    t0=MPI_Wtime()

#if defined (__icepack)
  call step_icepack(mesh, ice, time_evp, time_advec, time_therm) ! EVP, advection and thermodynamic parts    
#else     
    
    !___________________________________________________________________________
    ! Dynamics
    SELECT CASE (ice%whichEVP)
        CASE (0)
            if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'  
            call EVPdynamics(ice, partit, mesh)
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
    ! Advection part
    ! old FCT routines
    ! call ice_TG_rhs
    ! call ice_fct_solve
    ! call cut_off
    ! new FCT routines from Sergey Danilov 08.05.2018
#if defined (__oifs)
   do i=1,myDim_nod2D+eDim_nod2D
       ice_temp(i) = ice_temp(i)*a_ice(i)
   end do
#endif /* (__oifs) */

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_TG_rhs_div...'//achar(27)//'[0m'
    call ice_TG_rhs_div(ice, partit, mesh) 
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m' 
    call ice_fct_solve(ice, partit, mesh)
    
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
    call ice_update_for_div(ice, partit, mesh)
    
#if defined (__oifs)
    do i=1,myDim_nod2D+eDim_nod2D
        if (a_ice(i)>0.0_WP) ice_temp(i) = ice_temp(i)/a_ice(i)
    end do
#endif /* (__oifs) */

    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call cut_off...'//achar(27)//'[0m'
    call cut_off(ice, partit, mesh)
    
    if (use_cavity) call cavity_ice_clean_ma(ice, partit, mesh)
    t2=MPI_Wtime()
    
    !___________________________________________________________________________
    ! Thermodynamic part
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call thermodynamics...'//achar(27)//'[0m'
    call thermodynamics(ice, partit, mesh)
#endif /* (__icepack) */

    !___________________________________________________________________________
    do node=1,myDim_nod2D+eDim_nod2D
        if ( (ice%uice(node)/=0.0_WP .and. mesh%ulevels_nod2d(node)>1) .or. &
             (ice%vice(node)/=0.0_WP .and. mesh%ulevels_nod2d(node)>1) ) then
            write(*,*) " --> found cavity velocity /= 0.0_WP , ", mype
            write(*,*) " ulevels_nod2d(n) = ", mesh%ulevels_nod2d(node)
            write(*,*) " U_ice(n) = ", ice%uice(node)
            write(*,*) " V_ice(n) = ", ice%vice(node)
            write(*,*)
        end if 
    end do
    t3=MPI_Wtime()
    rtime_ice = rtime_ice + (t3-t0)
    rtime_tot = rtime_tot + (t3-t0)
    if(mod(step,logfile_outfreq)==0 .and. mype==0) then 
        write(*,*) '___ICE STEP EXECUTION TIMES____________________________'
#if defined (__icepack)
        write(*,"(A, ES10.3)") '    Ice Dyn.        :', time_evp
        write(*,"(A, ES10.3)") '        Ice Advect.     :', time_advec
        write(*,"(A, ES10.3)") '        Ice Thermodyn.  :', time_therm
#else
        write(*,"(A, ES10.3)") '    Ice Dyn.        :', t1-t0
        write(*,"(A, ES10.3)") '    Ice Advect.     :', t2-t1
        write(*,"(A, ES10.3)") '    Ice Thermodyn.  :', t3-t2
#endif /* (__icepack) */
        write(*,*) '   _______________________________'
        write(*,"(A, ES10.3)") '    Ice TOTAL       :', t3-t0
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
    use MOD_MESH
    use o_PARAM, only: WP   
    use g_CONFIG
    implicit none
    !
    type(t_ice)   , intent(inout), target :: ice
    type(t_tracer), intent(in)   , target :: tracers
    type(t_partit), intent(inout), target :: partit
    type(t_mesh)  , intent(in)   , target :: mesh
    !___________________________________________________________________________
    integer                               :: node
    !!PS character(MAX_PATH)                   :: filename
    !!PS real(kind=WP), external               :: TFrez  ! Sea water freeze temperature.
    !___________________________________________________________________________
    ! pointer on necessary derived types
    real(kind=WP), dimension(:)  , pointer:: a_ice, m_ice, m_snow, u_ice, v_ice
#include "associate_part_def.h"
#include "associate_mesh_def.h"
#include "associate_part_ass.h"
#include "associate_mesh_ass.h"
    a_ice  => ice%data(1)%values
    m_ice  => ice%data(2)%values
    m_snow => ice%data(3)%values
    u_ice  => ice%uice(:)
    v_ice  => ice%vice(:)
    
    !___________________________________________________________________________
    ! pointer on necessary derived types
    if(mype==0) write(*,*) 'initialize the sea ice'
    
    !___________________________________________________________________________
    do node=1,myDim_nod2D+eDim_nod2D    
    
        !_______________________________________________________________________
        if (ulevels_nod2d(node)>1) then
            !!PS m_ice(i)  = 1.0e15_WP
            !!PS m_snow(i) = 0.1e15_WP
            cycle ! --> if cavity, no sea ice, no initial state
        endif    
        
        !_______________________________________________________________________
        if (tracers%data(1)%values(1,node)< 0.0_WP) then
            if (geo_coord_nod2D(2,node)>0.0_WP) then
                m_ice(node) = 1.0_WP
                m_snow(node)= 0.1_WP 
            else
                m_ice(node) = 2.0_WP
                m_snow(node)= 0.5_WP 
            end if
            a_ice(node) = 0.9_WP
            u_ice(node) = 0.0_WP
            v_ice(node) = 0.0_WP
            
        endif
    enddo
end subroutine ice_initial_state
