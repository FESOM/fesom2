module ice_setup_step_interfaces
  interface
    subroutine ice_array_setup(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine

    subroutine ice_initial_state(mesh)
      use mod_mesh
      type(t_mesh), intent(in)  , target :: mesh
    end subroutine
  end interface
end module

!
!
!_______________________________________________________________________________
! ice initialization + array allocation + time stepping
subroutine ice_setup(mesh)
    use o_param
    use g_parsup
    use i_param
    use i_arrays
    use g_CONFIG
    use mod_mesh
    use ice_setup_step_interfaces
    implicit none 
    type(t_mesh), intent(in)           , target :: mesh
    
    ! ================ DO not change
    ice_dt=real(ice_ave_steps,WP)*dt
    ! ice_dt=dt
    Tevp_inv=3.0_WP/ice_dt 
    Clim_evp=Clim_evp*(evp_rheol_steps/ice_dt)**2/Tevp_inv  ! This is combination 
                                                            ! it always enters

    ! ================
    call ice_array_setup(mesh)
    call ice_fct_init(mesh)
    ! ================
    ! Initialization routine, user input is required 
    ! ================
    !call ice_init_fields_test
    call ice_initial_state(mesh)   ! Use it unless running test example
    if(mype==0) write(*,*) 'Ice is initialized'
end subroutine ice_setup
!
!
!_______________________________________________________________________________
subroutine ice_array_setup(mesh)
!
! inializing sea ice model 
!
! Variables that serve for exchange with atmosphere are nodal, to keep 
! back compatibility with FESOM input routines
 
use o_param
use i_param
use MOD_MESH
use i_arrays
use g_parsup
USE g_CONFIG

implicit none
type(t_mesh), intent(in)           , target :: mesh
integer   :: n_size, e_size, mn, k, n, n1, n2

#include  "associate_mesh.h"

n_size=myDim_nod2D+eDim_nod2D
e_size=myDim_elem2D+eDim_elem2D

! Allocate memory for variables of ice model
 allocate(u_ice(n_size), v_ice(n_size))
 allocate(U_rhs_ice(n_size), V_rhs_ice(n_size))
 allocate(sigma11(e_size), sigma12(e_size), sigma22(e_size))
 allocate(eps11(e_size),     eps12(e_size),   eps22(e_size))
 allocate(m_ice(n_size), a_ice(n_size), m_snow(n_size))
 allocate(rhs_m(n_size), rhs_a(n_size), rhs_ms(n_size))
 allocate(t_skin(n_size))
 allocate(U_ice_old(n_size), V_ice_old(n_size)) !PS
 allocate(m_ice_old(n_size), a_ice_old(n_size), m_snow_old(n_size), thdgr_old(n_size)) !PS
 if (whichEVP > 0) then
    allocate(u_ice_aux(n_size), v_ice_aux(n_size))
    allocate(alpha_evp_array(myDim_elem2D))
    allocate(beta_evp_array(n_size))

    alpha_evp_array=alpha_evp
    beta_evp_array =alpha_evp  ! alpha=beta works most reliable
    u_ice_aux=0.0_WP
    v_ice_aux=0.0_WP
 end if
 
 allocate(rhs_mdiv(n_size), rhs_adiv(n_size), rhs_msdiv(n_size))

 m_ice_old=0.0_WP !PS
 a_ice_old=0.0_WP !PS
 m_snow_old=0.0_WP !PS
 thdgr_old=0.0_WP !PS
 U_ice_old=0.0_WP !PS
 V_ice_old=0.0_WP !PS
 
 rhs_m=0.0_WP
 rhs_ms=0.0_WP
 rhs_a=0.0_WP
 m_ice=0.0_WP
 a_ice=0.0_WP
 m_snow=0.0_WP
 U_rhs_ice=0.0_WP
 V_rhs_ice=0.0_WP
 U_ice=0.0_WP
 V_ice=0.0_WP
 sigma11=0.0_WP
 sigma22=0.0_WP
 sigma12=0.0_WP
 eps11=0.0_WP
 eps12=0.0_WP
 eps22=0.0_WP
 t_skin=0.0_WP
 rhs_mdiv=0.0_WP
 rhs_adiv=0.0_WP
 rhs_msdiv=0.0_WP


! Allocate memory for arrays used in coupling 
! with ocean and atmosphere
 allocate(S_oc_array(n_size), T_oc_array(n_size))  ! copies of ocean T ans S
 S_oc_array = 0.0_WP
 T_oc_array = 0.0_WP
 allocate(fresh_wa_flux(n_size), net_heat_flux(n_size))
 fresh_wa_flux = 0.0_WP
 net_heat_flux = 0.0_WP
 allocate(stress_atmice_x(n_size), stress_atmice_y(n_size))    
 stress_atmice_x = 0.0_WP
 stress_atmice_y = 0.0_WP
 allocate(elevation(n_size))           ! =ssh  of ocean        
 elevation = 0.0_WP
 allocate(stress_iceoce_x(n_size), stress_iceoce_y(n_size))    
 stress_iceoce_x = 0.0_WP
 stress_iceoce_y = 0.0_WP
 allocate(U_w(n_size), V_w(n_size))   ! =uf and vf of ocean at surface nodes
#if defined (__oasis)
  allocate(oce_heat_flux(n_size), ice_heat_flux(n_size))
  allocate(tmp_oce_heat_flux(n_size), tmp_ice_heat_flux(n_size))
#if defined (__oifs)
  allocate(ice_alb(n_size), ice_temp(n_size), enthalpyoffuse(n_size))
  allocate(rhs_tempdiv(n_size), rhs_temp(n_size))
  ice_alb=0.6_WP
  ice_temp=265.15_WP
  rhs_tempdiv=0._WP
  rhs_temp=0._WP
  enthalpyoffuse=0._WP
#endif /* (__oifs) */
  oce_heat_flux=0._WP
  ice_heat_flux=0._WP
  tmp_oce_heat_flux=0._WP
  tmp_ice_heat_flux=0._WP
#endif /* (__oasis) */
end subroutine ice_array_setup
!
!
!
!_______________________________________________________________________________
! Sea ice model step
subroutine ice_timestep(step, mesh)
use i_arrays
use o_param
use g_parsup
use g_CONFIG
use i_PARAM, only: whichEVP
use mod_mesh

#if defined (__icepack)
    use icedrv_main,   only: step_icepack 
#endif

implicit none 
type(t_mesh), intent(in)   , target :: mesh
integer                    :: step,i
REAL(kind=WP)              :: t0,t1, t2, t3

#if defined (__icepack)
real(kind=WP)              :: time_evp, time_advec, time_therm
#endif

t0=MPI_Wtime()

#if defined (__icepack)
    call step_icepack(mesh, time_evp, time_advec, time_therm) ! EVP, advection and thermodynamic parts    
#else     
    
    !___________________________________________________________________________
    ! ===== Dynamics
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call EVPdynamics...'//achar(27)//'[0m'  
    SELECT CASE (whichEVP)
    CASE (0)
        call EVPdynamics(mesh)
    CASE (1)
        call EVPdynamics_m(mesh)
    CASE (2)
        call EVPdynamics_a(mesh)
    CASE DEFAULT
        if (mype==0) write(*,*) 'a non existing EVP scheme specified!'
        call par_ex
        stop
    END SELECT
    
    if (use_cavity) call cavity_ice_clean_vel(mesh)
    t1=MPI_Wtime()   
    
    !___________________________________________________________________________
    ! ===== Advection part
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
    call ice_TG_rhs_div(mesh)   
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_fct_solve...'//achar(27)//'[0m' 
    call ice_fct_solve(mesh)
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call ice_update_for_div...'//achar(27)//'[0m'
    call ice_update_for_div(mesh)
#if defined (__oifs)
    do i=1,myDim_nod2D+eDim_nod2D
        if (a_ice(i)>0.0_WP) ice_temp(i) = ice_temp(i)/a_ice(i)
    end do
#endif /* (__oifs) */
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call cut_off...'//achar(27)//'[0m'
    call cut_off(mesh)
    
    if (use_cavity) call cavity_ice_clean_ma(mesh)
    t2=MPI_Wtime()
    
    !___________________________________________________________________________
    ! ===== Thermodynamic part
    if (flag_debug .and. mype==0)  print *, achar(27)//'[36m'//'     --> call thermodynamics...'//achar(27)//'[0m'
    call thermodynamics(mesh)
#endif /* (__icepack) */


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
subroutine ice_initial_state(mesh)
    use i_ARRAYs
    use MOD_MESH
    use o_PARAM   
    use o_arrays        
    use g_parsup 
    use g_CONFIG
    implicit none
    !
    type(t_mesh), intent(in)           , target :: mesh
    integer                            :: i
    character(MAX_PATH)                      :: filename
    real(kind=WP), external            :: TFrez  ! Sea water freeze temperature.

#include  "associate_mesh.h"

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
        if (tr_arr(1,i,1)< 0.0_WP) then
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
