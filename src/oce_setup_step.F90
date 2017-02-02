! =================================================================
subroutine ocean_setup
USE o_PARAM
USE g_PARSUP
USE o_ARRAYS
USE g_config
IMPLICIT NONE
         
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
	end if    
        
	call stiff_mat
	if(mype==0) write(*,*) 'Arrays are set'
        
	!if(open_boundary) call set_open_boundary   !TODO
	
	if (tracer_adv==2) call fct_init
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
use g_config
use g_forcing_arrays
use o_mixing_kpp_mod ! KPP
IMPLICIT NONE
integer     :: elem_size, node_size

elem_size=myDim_elem2D+eDim_elem2D
node_size=myDim_nod2D+eDim_nod2D


! ================
! Velocities
! ================     
!allocate(stress_diag(2, elem_size))!delete me
allocate(UV(2, nl-1, elem_size))

if (use_means) allocate(UV_mean(2,nl-1, elem_size))
allocate(UV_rhs(2,nl-1, elem_size))
allocate(UV_rhsAB(2,nl-1, elem_size))
allocate(Visc(nl-1, elem_size))
! ================
! elevation and its rhs
! ================
allocate(eta_n(node_size), d_eta(node_size))
if (use_means) allocate(eta_n_mean(node_size))
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
if (use_means) allocate(Wvel_mean(nl, node_size))
! ================
! Temperature and salinity
! ================
allocate(T_rhs(nl-1, node_size))
allocate(S_rhs(nl-1, node_size))
allocate(tr_arr(nl-1,node_size,num_tracers),tr_arr_old(nl-1,node_size,num_tracers))
allocate(del_ttf(nl-1,node_size))
if (use_means) allocate(tr_arr_mean(nl-1,node_size,num_tracers))

allocate(bvfreq(nl,node_size),mixlay_dep(node_size),bv_ref(node_size))
! ================
! Ocean forcing arrays
! ================
allocate(Tclim(nl-1,node_size), Sclim(nl-1, node_size))
allocate(stress_surf(2,myDim_elem2D))    !!! Attention, it is shorter !!! 
allocate(relax2clim(node_size)) 
allocate(heat_flux(node_size), Tsurf(node_size))
allocate(water_flux(node_size), Ssurf(node_size))
! =================
! Arrays used to organize surface forcing
! =================
allocate(Tsurf_t(node_size,2), Ssurf_t(node_size,2))
allocate(tau_x_t(node_size,2), tau_y_t(node_size,2))  

! =================
! All auxiliary arrays
! =================
 
if(mom_adv==3) then
allocate(vorticity(nl-1,node_size))
vorticity=0.0_8
end if

! =================
! Visc and Diff coefs
! =================

allocate(Av(nl,elem_size), Kv(nl,node_size))
 
if (trim(mix_scheme)=='KPP') then
   allocate(Kv2(nl,node_size,num_tracers))
   Kv2=0d0 
else if(trim(mix_scheme)=='PP') then   
   if (AvKv) then
      allocate(Kv2(nl,node_size,2))
      Kv2=0d0 
   end if 
else
   stop("!not existing mixing scheme!")
   call par_ex
end if

if (use_means .and. hbl_diag) then
   allocate(hbl_mean(node_size)) 
end if

if (use_means .and. AvKv) then
   allocate(Av_mean(nl,elem_size),Kv_mean(nl,node_size,num_tracers)) 
end if

Av=0.0_WP
Kv=0.0_WP

if (trim(mix_scheme)=='KPP') call oce_mixing_kpp_init ! Setup constants, allocate arrays and construct look up table

!Velocities at nodes
allocate(Unode(2,nl-1,node_size))

!Tracer gradients&RHS  
allocate(ttrhs(nl-1,node_size))
allocate(tr_xy(2,nl-1,myDim_elem2D+eDim_elem2D+eXDim_elem2D))
allocate(neutral_slope(3, nl-1, node_size))
allocate(Kd(4, nl-1, myDim_nod2D+eDim_nod2D))
neutral_slope=0.0_WP

allocate(sigma_xy(2, nl-1, node_size))
sigma_xy=0.0_WP
Kd=0.0_WP

allocate(sw_beta(nl-1, node_size), sw_alpha(nl-1, node_size))
sw_beta=0.0_WP
sw_alpha=0.0_WP

if (Fer_GM) then
   allocate(fer_c(node_size), fer_gamma(2, nl, node_size), fer_K(node_size))
   allocate(fer_wvel(nl, node_size), fer_UV(2, nl-1, elem_size))
   if (use_means) then
      allocate(fer_UV_mean(2, nl-1, elem_size), fer_wvel_mean(nl, node_size))
   end if
!   allocate(sw_beta(nl-1, node_size), sw_alpha(nl-1, node_size))  KPP needs this even if Fer_GM is not called
!   sw_beta=0.0_WP
!   sw_alpha=0.0_WP
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

    S_rhs=0.0_WP
    water_flux=0.0_WP
    Ssurf=0.0_WP

    tr_arr=0d0
    tr_arr_old=0d0    

    bvfreq=0d0
    mixlay_dep=0d0
    bv_ref=0d0

    relax2clim=0.0_WP

    Tsurf_t=0.0_WP
    Ssurf_t=0.0_WP
    tau_x_t=0.0_WP
    tau_y_t=0.0_WP

    if (use_means) then
       UV_mean=0.0_WP
       Wvel_mean=0.0_WP
       eta_n_mean=0.0_WP
       tr_arr_mean=0.0_WP
       if (Fer_GM) then
          fer_UV_mean=0.0_WP
          fer_wvel_mean=0.0_WP
       end if
       if(AvKv) then
          Av_mean=0d0
          Kv_mean=0d0
       endif
       if (hbl_diag) then
          hbl_mean=0d0
       endif
    endif
END SUBROUTINE array_setup
!==========================================================================
SUBROUTINE oce_timestep(n)
use g_config, only: logfile_outfreq
USE o_MESH
USE o_ARRAYS
USE o_PARAM
USE g_PARSUP
use o_mixing_kpp_mod ! _OG_


use g_comm_auto !??????????????????????????????
use ieee_arithmetic        !??????????????????????????????


IMPLICIT NONE
real(kind=8)      :: t1, t2, t3, t4, t5

integer i,elem,nz,n
real*8,external:: dnrm2 
real*8         :: global_max_hf ,global_min_hf
real(kind=8), allocatable :: u_aux3(:,:), v_aux3(:,:), aux1(:)
real(kind=8)      :: global_vol_eta ,local_vol_eta ,global_max_eta ,global_min_eta ,local_max_eta ,local_min_eta
integer tr_num
integer           :: el !??????????????????????????????
CHARACTER(LEN=80) :: fmt   !??????????????????????????????


  t1=MPI_Wtime() 
  call pressure_bv

  call status_check



if (trim(mix_scheme)=='KPP') then
  call oce_mixing_KPP(Av, Kv2)
! oce_mixing_kpp.F90 should be modified in a way that works for each tracer for the case double diffusion
! tr_num = 1 temp
! tr_num = 2 salt 
!   if (tr_num==1) then 
!      Kv2(:,:,1)=Kv2(:,:,2)
      Kv(:,:)=Kv2(:,:,1)

else if(trim(mix_scheme)=='PP') then
  call oce_mixing_PP
  if (AvKv) then
     Kv2(:,:,1)=Kv(:,:)  ! for the output 
     Kv2(:,:,2)=Kv(:,:)  ! for the output
  end if 
else
   stop("!not existing mixing scheme!")
   call par_ex
end if  

  if(mom_adv/=3) then
    call compute_vel_rhs
  else
    call compute_vel_rhs_vinv
  end if
  if (tau_c > 1.e-12) call viscosity_filt2x

  if (i_vert_visc)    call impl_vert_visc ! Probably should be moved for Btr-bcl splitting case

  call compute_ssh_rhs
  t2=MPI_Wtime()
  call solve_ssh
  t3=MPI_Wtime() 

  call update_vel

  if (Fer_GM) then
     call fer_compute_C_K
     call compute_sigma_xy(tr_arr(:,:,1),tr_arr(:,:,2))
     call fer_solve_Gamma
     call fer_gamma2vel
  end if

  call vert_vel
  t4=MPI_Wtime()
  call solve_tracers
  t5=MPI_Wtime() 
  if(mype==0) write(*,*) 'SSH PE=0 ', maxval(eta_n), minval(eta_n)

  if(mype==0) then  
   write(*,*) 'Step took   ', t5-t1
   write(*,*) 'Solver      ', t3-t2
   write(*,*) 'Dynamics    ', t2-t1
   write(*,*) 'Update+W    ', t4-t3
   write(*,*) 'Tracer      ', t5-t4
  end if
  !___________________________________________________________________________
	! write out field estimates
	if(mod(n,logfile_outfreq)==0 .and. mype==0) then
! 		fmt="(A, ES10.3, A, ES10.3)"
		fmt="(A, F10.3, A, F10.3)"
		write(*,*) '	___ALE OCEAN STEP: FIELDS______________________________'
		write(*,fmt) '	min(eta_n) = ',minval(eta_n)          ,', max(eta_n) = ', maxval(eta_n)
		write(*,fmt) '	min(U)     = ',minval(UV(1,:,:))      ,', max(U)     = ', maxval(UV(1,:,:))
		write(*,fmt) '	min(V)     = ',minval(UV(2,:,:))      ,', max(V)     = ', maxval(UV(2,:,:))
		write(*,fmt) '	min(Wvel)  = ',minval(Wvel(:,:))      ,', max(Wvel)  = ', maxval(Wvel(:,:))
		write(*,fmt) '	min(temp)  = ',minval(tr_arr(:,:,1))  ,', max(temp)  = ', maxval(tr_arr(:,:,1))
		write(*,fmt) '	min(salt)  = ',minval(tr_arr(:,:,2))  ,', max(salt)  = ', maxval(tr_arr(:,:,2))
		write(*,*)
		write(*,fmt) '	min(hflux) = ',minval(heat_flux)      ,'  max(hflux) = ', maxval(heat_flux)
		write(*,fmt) '	min(wflux) = ',minval(water_flux)     ,'  max(wflux) = ', maxval(water_flux)
		write(*,*)
	end if
	
	
	!___________________________________________________________________________
	! check integral of global ssh
	if (mod(n,logfile_outfreq)==0) then
		global_max_eta=0.0_WP
		global_min_eta=0.0_WP
		local_max_eta=maxval(heat_flux)
		local_min_eta=minval(heat_flux)
		
		call MPI_AllREDUCE(local_max_eta, global_max_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
							MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(local_min_eta, global_min_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
							MPI_COMM_WORLD, MPIerr)
							
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,"(A, I, A, ES10.3, A, ES10.3)") ' --mstep=',mstep,'--> global min/max hflux      : ',global_min_eta,' / ',global_max_eta		
			write(*,*)
		endif
	
		local_vol_eta  = local_vol_eta/myDim_elem2D
		local_max_eta  = maxval(eta_n)
		local_min_eta  = minval(eta_n)
		global_vol_eta=0.0_WP
		global_max_eta=0.0_WP
		global_min_eta=0.0_WP
		call MPI_AllREDUCE(local_vol_eta, global_vol_eta, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
							MPI_COMM_WORLD, MPIerr)
		global_vol_eta = global_vol_eta/npes
		call MPI_AllREDUCE(local_max_eta, global_max_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
							MPI_COMM_WORLD, MPIerr)
		call MPI_AllREDUCE(local_min_eta, global_min_eta, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
							MPI_COMM_WORLD, MPIerr)
		
		if (mod(n,logfile_outfreq)==0 .and. mype==0) then
			write(*,*)
			write(*,"(A, I, A, ES10.3, A, ES10.3, A, ES10.3)") ' --mstep=',n,'--> global max/min/int eta_n : ',global_max_eta,' / ',global_min_eta,' / ',global_vol_eta
		endif
	endif 
		
	
! 	call gather_nod(Wvel(1,:),aux1)
! 	if (mype==0) then
! 		write(*,*) ' --CHECK--> global mean Wvel: ',sum(aux1)/nod2d
! 		write(*,*) ' --CHECK--> global max  Wvel: ',maxval(aux1)
! 		write(*,*) ' --CHECK--> global mean Wvel: ',minval(aux1)
! 		write(*,*)
! 	endif
	!___________________________________________________________________________
	! check tracerfor nan
	do el=1,myDim_nod2d
		do nz=1,nl
			if ( ieee_is_nan(tr_arr(nz, el,1)) .or. &
			     tr_arr(nz, el,1) < -5.0 .or. tr_arr(nz, el,1)>100) then
				write(*,*) ' --STOP--> found temperture become NaN or <-5.0, >100'
				write(*,*) 'mype        = ',mype
				write(*,*) 'n           = ',n
				write(*,*) 'el          = ',el
				write(*,*) 'nz          = ',nz
				write(*,*) 'temp(nz, el)) = ',tr_arr(nz, el,1)
				write(*,*) 'hflux       = ',heat_flux(el)
		
				call par_ex(1)
			endif
			if ( ieee_is_nan(tr_arr(nz, el,2)) .or.  &
			     tr_arr(nz, el,2) < 0 .or. tr_arr(nz, el,2)>100 ) then
				write(*,*) ' --STOP--> found salinity become NaN or <0, >100'
				write(*,*) 'mype        = ',mype
				write(*,*) 'n           = ',n
				write(*,*) 'el          = ',el
				write(*,*) 'nz          = ',nz
				write(*,*) 'salt(nz, el)) = ',tr_arr(nz, el,2)
				write(*,*) 'wflux       = ',water_flux(el)
				call par_ex(1)
			endif 
		end do
	end do
END SUBROUTINE oce_timestep
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
  
  call read_init_ts
     Tclim=tr_arr(:,:,1)
     Sclim=tr_arr(:,:,2)
     Tsurf=tr_arr(1,:,1)
     Ssurf=tr_arr(1,:,2)
  if(mype==0) write(*,*) 'read T/S climatology', trim(OceClimaDataName)
  !=====================
  ! Restart from NetCDF output (not fully rigorous because of AB)
  !=====================    
  if (r_restart) then
     if(mype==0) write(*,*) 'read ocean restart file'
     call oce_input
  end if

  relax2clim=0.0 
end subroutine oce_initial_state
!==========================================================================
