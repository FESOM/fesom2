!
! ice initialization + array allocation + time stepping
!
!==============================================================================
subroutine ice_setup
use o_param
use g_parsup
use o_mesh
use i_param
use i_arrays
use g_CONFIG
implicit none 
  ! ================ DO not change
  ice_dt=ice_ave_steps*dt
 ! ice_dt=dt
  Tevp_inv=3.0/ice_dt 
  Clim_evp=Clim_evp*(evp_rheol_steps/ice_dt)**2/Tevp_inv  ! This is combination 
                                                      ! it always enters 
  ! ================
  call ice_array_setup
  call ice_fct_init 
  ! ================
  ! Initialization routine, user input is required 
  ! ================
  !call ice_init_fields_test
  call ice_initial_state     ! Use it unless running test example
  if(mype==0) write(*,*) 'Ice is initialized'
end subroutine ice_setup
! ============================================================================
subroutine ice_array_setup
!
! inializing sea ice model 
!
! Variables that serve for exchange with atmosphere are nodal, to keep 
! back compatibility with FESOM input routines
 

use o_param
use i_param
use o_mesh
use i_arrays
use g_parsup
USE g_CONFIG

implicit none
integer   :: n_size, e_size, mn, k, n, n1, n2

n_size=myDim_nod2D+eDim_nod2D
e_size=myDim_elem2D+eDim_elem2D

! Allocate memory for variables of ice model      
 allocate(u_ice(n_size), v_ice(n_size))
 if (use_means) allocate(u_ice_mean(n_size), v_ice_mean(n_size))
 allocate(U_rhs_ice(n_size), V_rhs_ice(n_size))
 allocate(sigma11(e_size), sigma12(e_size), sigma22(e_size)) 
 allocate(m_ice(n_size), a_ice(n_size), m_snow(n_size))
 if (use_means) allocate(m_ice_mean(n_size), a_ice_mean(n_size), m_snow_mean(n_size))
 allocate(rhs_m(n_size), rhs_a(n_size), rhs_ms(n_size))
 allocate(t_skin(n_size))

 allocate(U_ice_old(n_size), V_ice_old(n_size)) !PS
 allocate(m_ice_old(n_size), a_ice_old(n_size), m_snow_old(n_size), thdgr_old(n_size)) !PS
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
 t_skin=0.0_WP

if (use_means) then
 m_ice_mean=0.0_WP
 a_ice_mean=0.0_WP
 m_snow_mean=0.0_WP
 U_ice_mean=0.0_WP
 V_ice_mean=0.0_WP
endif

! Allocate memory for arrays used in coupling 
! with ocean and atmosphere
 allocate(S_oc_array(n_size), T_oc_array(n_size))  ! copies of ocean T ans S
 allocate(fresh_wa_flux(n_size), net_heat_flux(n_size))
 allocate(stress_atmice_x(n_size), stress_atmice_y(n_size))    
 allocate(stress_atmoce_x(n_size), stress_atmoce_y(n_size))    
 allocate(elevation(n_size))           ! =ssh  of ocean        
 allocate(stress_iceoce_x(n_size), stress_iceoce_y(n_size))    
 allocate(U_w(n_size), V_w(n_size))   ! =uf and vf of ocean at surface nodes
end subroutine ice_array_setup
!==========================================================================
subroutine ice_timestep(step)
! 
! Sea ice model step
!
use o_param
use g_parsup
USE g_CONFIG
implicit none
integer      :: step 
REAL(kind=8) :: t0,t1, t2, t3
t0=MPI_Wtime()
 ! ===== Dynamics
 call EVPdynamics
 t1=MPI_Wtime()     
 ! ===== Advection part
 call ice_fct_solve
 call cut_off
 t2=MPI_Wtime()
 ! ===== Thermodynamic part
 call thermodynamics
 t3=MPI_Wtime()
if(mod(step,logfile_outfreq)==0 .and. mype==0) then 
		write(*,*) '___ICE STEP EXECUTION TIMES____________________________'
		write(*,"(A, ES10.3)") '	Ice Dyn.        :', t1-t0
		write(*,"(A, ES10.3)") '	Ice Advect.     :', t2-t1
		write(*,"(A, ES10.3)") '	Ice Thermodyn.  :', t3-t2
		write(*,*) '   _______________________________'
		write(*,"(A, ES10.3)") '	Ice TOTAL       :', t3-t0
		write(*,*)
endif

end subroutine ice_timestep
!==============================================================================

subroutine ice_initial_state
  !sets inital values or reads restart file for ice model
  use i_ARRAYs
  use o_MESH    
  use o_PARAM   
  use o_arrays        
  use g_parsup 
  USE g_CONFIG
  implicit none
  !
  integer        :: i
  character*100  :: filename
  real*8, external  :: TFrez  ! Sea water freeze temperature.
  m_ice =0.
  a_ice =0.
  u_ice =0.
  v_ice =0.
  m_snow=0.
  if(mype==0) write(*,*) 'initialize the sea ice'

  do i=1,myDim_nod2D+eDim_nod2D                           
     if (tr_arr(1,i,1)< 0.0_WP) then
	if (geo_coord_nod2D(2,i)>0.) then
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
