! =================================================================
subroutine ocean_setup
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_config
USE g_forcing_param, only: use_virt_salt
IMPLICIT NONE

       if (use_ALE) then
	  !___setup virt_salt_flux____________________________________________________
          ! if the ale thinkness remain unchanged (like in 'linfs' case) the vitrual 
          ! salinity flux need to be used
          ! otherwise we set the reference salinity to zero
          if ( .not. trim(which_ALE)=='linfs') then
 	     use_virt_salt=.false.
             ! this will force the virtual saltinity flux to be zero
             ref_sss_local=.false.
             ref_sss=0._WP
             is_nonlinfs = 1.0_WP
          else
             use_virt_salt=.true.
             is_nonlinfs = 0.0_WP
          end if
        end if
        call array_setup
    !___________________________________________________________________________
	! initialize arrays for ALE
	if(use_ALE) then
		if(mype==0) then
				write(*,*) '____________________________________________________________'
				write(*,*) ' --> call ale_init'
				write(*,*)
		end if	
		call ale_init
		call stiff_mat_ale
	else
		call stiff_mat
	end if    
        
	
	if(mype==0) write(*,*) 'Arrays are set'
        
	!if(open_boundary) call set_open_boundary   !TODO
	
	call fct_init
        call muscl_adv_init
	!=====================
	! Initialize fields
	! A user-defined routine has to be called here!
	!=====================
	if(toy_ocean) then  
#ifdef NA_TEST
	 call init_fields_na_test
#else
	 call initial_state_test
#endif
	 !call initial_state_channel_test 
	 !call initial_state_channel_narrow_test
	 !call init_fields_na_test  
 	 !call init_fields_global_test
        else
#ifdef NA_TEST
	 call init_fields_na_test
#else
         call oce_initial_state   ! Use it if not running tests
#endif
        end if

         if (.not.r_restart) tr_arr_old=tr_arr
         
    !___________________________________________________________________________
	! first time fill up array for hnode & helem
	if(use_ALE) then
		if(mype==0) then
				write(*,*) '____________________________________________________________'
				write(*,*) ' --> call init_thickness_ale'
				write(*,*)
		end if	
		call init_thickness_ale
	end if     

	 if(mype==0) write(*,*) 'Initial state'
if (w_split .and. mype==0) then
	write(*,*) '******************************************************************************'
	write(*,*) 'vertical velocity will be split onto explicit and implicit constitutes;'
	write(*,*) 'maximum explicit W is set to: ', w_exp_max
	write(*,*) '******************************************************************************'
end if
end subroutine ocean_setup

!==========================================================
!
SUBROUTINE array_setup
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use g_comm_auto
use g_config
use g_forcing_arrays
use o_mixing_kpp_mod ! KPP
USE g_forcing_param, only: use_virt_salt
IMPLICIT NONE
integer     :: elem_size, node_size
integer     :: n
elem_size=myDim_elem2D+eDim_elem2D
node_size=myDim_nod2D+eDim_nod2D


! ================
! Velocities
! ================     
!allocate(stress_diag(2, elem_size))!delete me
allocate(UV(2, nl-1, elem_size))
allocate(UV_rhs(2,nl-1, elem_size))
allocate(UV_rhsAB(2,nl-1, elem_size))
allocate(Visc(nl-1, elem_size))
! ================
! elevation and its rhs
! ================
allocate(eta_n(node_size), d_eta(node_size))
allocate(ssh_rhs(node_size))
! ================
! Monin-Obukhov
! ================
if (use_ice .and. mo_on) allocate(mo(nl,node_size),mixlength(node_size))
if (use_ice .and. mo_on) mixlength=0.
! ================
! Vertical velocity and pressure
! ================
allocate(Wvel(nl, node_size), hpressure(nl,node_size))
allocate(Wvel_e(nl, node_size), Wvel_i(nl, node_size))
! ================
! Temperature and salinity
! ================
allocate(T_rhs(nl-1, node_size))
allocate(S_rhs(nl-1, node_size))
allocate(tr_arr(nl-1,node_size,num_tracers),tr_arr_old(nl-1,node_size,num_tracers))
allocate(del_ttf(nl-1,node_size))
allocate(bvfreq(nl,node_size),mixlay_dep(node_size),bv_ref(node_size))
! ================
! Ocean forcing arrays
! ================
allocate(Tclim(nl-1,node_size), Sclim(nl-1, node_size))
allocate(stress_surf(2,myDim_elem2D))    !!! Attention, it is shorter !!! 
allocate(relax2clim(node_size)) 
allocate(heat_flux(node_size), Tsurf(node_size))
allocate(water_flux(node_size), Ssurf(node_size))
allocate(relax_salt(node_size))
allocate(virtual_salt(node_size))

allocate(heat_flux_old(node_size), Tsurf_old(node_size)) !PS
allocate(water_flux_old(node_size), Ssurf_old(node_size)) !PS
allocate(real_salt_flux(node_size)) !PS
! =================
! Arrays used to organize surface forcing
! =================
allocate(Tsurf_t(node_size,2), Ssurf_t(node_size,2))
allocate(tau_x_t(node_size,2), tau_y_t(node_size,2))  

! =================
! All auxiliary arrays
! =================
 
!if(mom_adv==3) then
allocate(vorticity(nl-1,node_size))
vorticity=0.0_8
!end if

! =================
! Visc and Diff coefs
! =================

allocate(Av(nl,elem_size), Kv(nl,node_size))

Av=0.0_WP
Kv=0.0_WP
if (trim(mix_scheme)=='KPP') then
   allocate(Kv_double(nl,node_size,num_tracers))
   Kv_double=0.0_WP
   call oce_mixing_kpp_init ! Setup constants, allocate arrays and construct look up table
end if

!Velocities at nodes
allocate(Unode(2,nl-1,node_size))

! tracer gradients & RHS  
allocate(ttrhs(nl-1,node_size))
allocate(tr_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(tr_z(nl,myDim_nod2D+eDim_nod2D))

! neutral slope etc. to be used in Redi formulation
allocate(neutral_slope(3, nl-1, node_size))
allocate(slope_tapered(3, nl-1, node_size))
allocate(Ki(node_size))
neutral_slope=0.0_WP
slope_tapered=0.0_WP

allocate(MLD1(node_size), MLD2(node_size), MLD1_ind(node_size), MLD2_ind(node_size))

do n=1, node_size
!  Ki(n)=K_hor*area(1,n)/scale_area
   Ki(n)=K_hor*mesh_resolution(n)/100000.0_WP
end do
call exchange_nod(Ki)

! xy gradient of a neutral surface
allocate(sigma_xy(2, nl-1, node_size))
sigma_xy=0.0_WP

! alpha and beta in the EoS
allocate(sw_beta(nl-1, node_size), sw_alpha(nl-1, node_size))
sw_beta=0.0_WP
sw_alpha=0.0_WP

if (Fer_GM) then
   allocate(fer_c(node_size), fer_gamma(2, nl, node_size), fer_K(node_size))
   allocate(fer_wvel(nl, node_size), fer_UV(2, nl-1, elem_size))
   fer_gamma=0.0_WP
   fer_uv=0.0_WP
   fer_wvel=0.0_WP
   fer_K=500.
   fer_c=1.
end if

! =================
! Initialize with zeros 
! =================

    UV=0.0_WP
    UV_rhs=0.0_WP
    UV_rhsAB=0.0_WP
!
    eta_n=0.0_WP
    ssh_rhs=0.0_WP
    Wvel=0.0_WP
    Wvel_e	=0.0_WP
    Wvel_i	=0.0_WP
    hpressure=0.0_WP
!
    T_rhs=0.0_WP
    heat_flux=0.0_WP
    Tsurf=0.0_WP
    heat_flux_old=0.0_WP !PS
    Tsurf_old=0.0_WP !PS

    S_rhs=0.0_WP
    water_flux=0.0_WP
    relax_salt=0.
    virtual_salt=0.

    Ssurf=0.0_WP
    water_flux_old=0.0_WP !PS
    Ssurf_old=0.0_WP !PS
    
    real_salt_flux=0.0_WP
    
    tr_arr=0d0
    tr_arr_old=0d0    

    bvfreq=0d0
    mixlay_dep=0d0
    bv_ref=0d0

    MLD1   =0.0_WP
    MLD2   =0.0_WP
    MLD1_ind=0
    MLD2_ind=0

    relax2clim=0.0_WP

    Tsurf_t=0.0_WP
    Ssurf_t=0.0_WP
    tau_x_t=0.0_WP
    tau_y_t=0.0_WP
END SUBROUTINE array_setup
!==========================================================================
SUBROUTINE oce_initial_state
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_config
USE g_input
  !
  ! reads the initial state or the restart file for the ocean
  !
  implicit none
  integer       :: i, node

  ! ===============
  ! read ocean state
  ! ===============
   if(mype==0) write(*,*) 'read T/S climatology', trim(OceClimaDataName)
   call read_init_ts
!  call ini_global_ocean ! initialize T&S somehow differently for a particular test case !
  
     Tclim=tr_arr(:,:,1)
     Sclim=tr_arr(:,:,2)
     Tsurf=tr_arr(1,:,1)
     Ssurf=tr_arr(1,:,2)
  relax2clim=0.0 
end subroutine oce_initial_state
!==========================================================================
