!=======================================================================
subroutine oce_mixing_pp
    !  Compute Richardson number dependent Av and Kv following
    !  Pacanowski and Philander, 1981
    !  Av = Avmax * factor**2 + Av0,
    !  Kv = Kvmax * factor**3 + Kv0, 
    !  factor=1/(1+5Ri),  Ri=N**2/(dU/dz)**2 is the Richardson number
    !                     N is the buoyancy frequency
    !                     dU/dz is the vertical velocity shear
    !  Avmax, Kvmax are tunable
    
    ! Output: Av(2:nlevels(elem)-1,:)  == vert. visc. coeff. at zbar(2:...)
    !         Kv(2:nlevels_nod2D(node),:) == vert. diff. coeff. at zbar(2:...)   
    ! NR external if for mo
    ! SD no if in Kv computations (only minor differences are introduced)
    !    
    !      
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_config
use i_arrays
IMPLICIT NONE

real(kind=WP)         :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s, Kv0_b
integer               :: node, nz, nzmax, elem, elnodes(3), i
real(kind=WP)         :: rhopot(nl), bulk_0(nl), bulk_pz(nl), bulk_pz2(nl)
real(kind=WP)         :: bulk_up, bulk_dn
real(kind=WP)         :: wndmix=1.e-3, wndnl=2, kv_conv=0.1_WP, av_conv=0.1_WP

	!___________________________________________________________________________
	if(use_ALE) then
		do node=1, myDim_nod2D+eDim_nod2D
			!___________________________________________________________________
			! implement changing ALE zlevel at every node in PP mxing 
			do nz=2,nlevels_nod2d(node)-1
				dz_inv=1.0_WP/(Z_3d_n(nz-1,node)-Z_3d_n(nz,node))
				shear = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
						(Unode(2,nz-1,node)-Unode(2,nz,node))**2 
				shear = shear*dz_inv*dz_inv
				Kv(nz,node) = shear/(shear+5.*max(bvfreq(nz,node),0.0_WP)+1.0e-14)  ! To avoid NaNs at start
			end do
		end do
		!_______________________________________________________________________
		! PP mixing withour ALE
	else
		do node=1, myDim_nod2D+eDim_nod2D
			!___________________________________________________________________
			do nz=2,nlevels_nod2d(node)-1
				dz_inv=1.0_WP/(Z(nz-1)-Z(nz))
				shear = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
						(Unode(2,nz-1,node)-Unode(2,nz,node))**2 
				shear = shear*dz_inv*dz_inv
				Kv(nz,node) = shear/(shear+5.*max(bvfreq(nz,node),0.0_WP)+1.0e-14)  ! To avoid NaNs at start
			end do
		enddo
	endif
	
	!___________________________________________________________________________
	! add vertical mixing scheme of Timmermann and Beckmann, 2004,"Parameterization 
	! of vertical mixing in the Weddell Sea!
	! Computes the mixing length derived from the Monin-Obukhov length
	! --> in FESOM1.4 refered as TB04 mixing scheme
	if (use_ice .and. mo_on) then !stress is only partial!!
		!_______________________________________________________________________
		if(use_ALE) then
			do node=1, myDim_nod2D+eDim_nod2D
				!_______________________________________________________________
				! calcualte monin obukov length
				call mo_length(water_flux(node),heat_flux(node), &         
								stress_atmoce_x(node),stress_atmoce_y(node), &    
								u_ice(node),v_ice(node),a_ice(node), &                             
								dt, mixlength(node))
				!_______________________________________________________________
				! increase vertical diffusion within monin obukov length to namelist
				! parameter modiff. modiff in moment set to 0.01 --> that means 
				! very strong vertical mixing within mixlength
				do nz = 2,nlevels_nod2d(node)-1
					mo(nz,node) = 0._WP
					if(abs(zbar_3d_n(nz,node)) <= mixlength(node)) mo(nz,node)=modiff    ! Potentialy bad place 
				end do 
			end do
		!_______________________________________________________________________
		else
			do node=1, myDim_nod2D+eDim_nod2D
				!_______________________________________________________________
				! calcualte monin obukov length
				call mo_length(water_flux(node),heat_flux(node), &         
								stress_atmoce_x(node),stress_atmoce_y(node), &    
								u_ice(node),v_ice(node),a_ice(node), &                             
								dt, mixlength(node))
				!_______________________________________________________________
				! increase vertical diffusion within monin obukov length to namelist
				! parameter modiff
				do nz = 2,nlevels_nod2d(node)-1
					mo(nz,node) = 0._WP
					if(abs(zbar(nz)) <= mixlength(node)) mo(nz,node)=modiff    ! Potentialy bad place 
				end do                                                        ! IF inside the internal cycle
			end do
		endif
		
		!_______________________________________________________________________
		! viscosity
		DO elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:,elem)
			DO nz=2,nlevels(elem)-1
				Av(nz,elem)= mix_coeff_PP*sum(Kv(nz,elnodes)**2)/3.0_WP + A_ver  &
							+ sum(mo(nz,elnodes))/3.0 
			END DO
		END DO
		!_______________________________________________________________________
		! diffusivity
		do node=1,myDim_nod2D+eDim_nod2D
			!___________________________________________________________________
			! set constant background diffusivity with namelist value K_ver
			if(Kv0_const) then
				do nz=2,nlevels_nod2d(node)-1 
					Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3+K_ver + mo(nz,node)
				end do
			!___________________________________________________________________________
			! set latitudinal and depth dependent background diffusivity after 
			! Q. Wang FESOM1.4 approach
			else
				do nz=2,nlevels_nod2d(node)-1 
					call Kv0_background_qiang(Kv0_b,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node)))
					Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3+Kv0_b + mo(nz,node)
				end do
			end if
		end do
	!___________________________________________________________________________
	! Original Pacanowski and Philander (PP) mixing scheme
	else ! .not.  (use_ice .and. mo_on)
		!_______________________________________________________________________
		! viscosity
		DO elem=1, myDim_elem2D
			elnodes=elem2D_nodes(:,elem)
			DO nz=2,nlevels(elem)-1
				Av(nz,elem)= mix_coeff_PP*sum(Kv(nz,elnodes)**2)/3.0_WP+A_ver
			END DO
		END DO
		!_______________________________________________________________________
		! diffusivity
		do node=1,myDim_nod2D+eDim_nod2D
			!___________________________________________________________________
			! set constant background diffusivity with namelist value K_ver
			if(Kv0_const) then
				do nz=2,nlevels_nod2d(node)-1 
					Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3+K_ver
				end do
			!___________________________________________________________________________
			! set latitudinal and depth dependent background diffusivity after 
			! Q. Wang FESOM1.4 approach
			else
				do nz=2,nlevels_nod2d(node)-1 
					call Kv0_background_qiang(Kv0_b,geo_coord_nod2D(2,node)/rad,abs(zbar_3d_n(nz,node)))
					Kv(nz,node) = mix_coeff_PP*Kv(nz,node)**3+Kv0_b
				end do
			end if
		end do
	end if
	
	!___________________________________________________________________________
	DO node=1, myDim_nod2D+eDim_nod2D
		DO nz=2, nlevels_nod2d(node)-1
			if (bvfreq(nz,node) < 0.) Kv(nz,node)=kv_conv
			if (nz<=wndnl+1) then
				Kv(nz,node)=max(Kv(nz,node), wndmix)
			end if
		END DO
	END DO
	
	DO elem=1, myDim_elem2D
		elnodes=elem2D_nodes(:,elem)
		DO nz=2,nlevels(elem)-1
			if (any(bvfreq(nz, elnodes) < 0.)) Av(nz,elem)=av_conv
			if (nz<=wndnl+1) then
				Av(nz,elem)=max(Av(nz,elem), wndmix)
			end if
		END DO
	END DO
	
end subroutine oce_mixing_pp
! ========================================================================
subroutine mo_length(water_flux,heat_flux,stress_x,stress_y,  &
       u_ice,v_ice,a_ice,dt,mixlength)
    ! vertical mixing scheme of Timmermann and Beckmann, 2004.
    ! computes the mixing length derived from the Monin-Obukhov length
    ! Ralph Timmermann, 14.06.2006
  
  USE o_PARAM, only: WP
    implicit none

    real(kind=WP)             :: water_flux, heat_flux, stress_x, stress_y
    real(kind=WP)             :: u_ice, v_ice, a_ice, uabs
    real(kind=WP)             :: dt, ret, rtc, mixlength
    real(kind=WP)             :: qfm, qtm, qw
    real(kind=WP)             :: ustar,tau, obuk
    real(kind=WP), parameter :: cosgam = 0.913632  ! cos(24.*3.14/180.)

    qfm            = water_flux * 34.       ! note that water_flux>0
    ![psu * m/s]   [m/s]   [psu]    ! f. upward fresh water flux

    qtm            = - 2.38e-7 * heat_flux  ! heat_flux>0 f. upward heat flux
    ![K * m/s]

    tau = sqrt(stress_x**2+stress_y**2)
    ustar = sqrt(tau/1030.)
    uabs = sqrt(u_ice**2+v_ice**2)

    qw = 1.25 * ustar**3 * (1.-a_ice) + 0.005 * uabs**3 * cosgam * a_ice  !Eq. 8 of TB04

    call pmlktmo(qfm,qtm,qw,obuk)

    rtc=dt/(10.0 * 86400.0)     ! (NR: inverse of) time constant of mixed layer retreat

    if (obuk.lt.mixlength) then
       ret = (obuk-mixlength)*rtc
       mixlength = mixlength+ret
    else
       mixlength=obuk
    endif

end subroutine mo_length
  !
  !=========================================================================	    	    
  !
subroutine pmlktmo(qfm,qtm,qw,obuk)
    ! gives the Monin-Obukhov length
    ! qtm  = Heat Flux into ML                                 [K m/s]
    ! qfm  = salinity flux into ML                             [psu m/s]
    ! qw   = production of turbulent kinetic energy 
    !-----------------------------------------------------------------------   
  USE o_PARAM, only: WP
    implicit none

    integer           :: iter
    real(kind=WP)            :: qtm, qfm, qw, obuk
    real(kind=WP), parameter :: qhw   = 1/7.0                 ! [1/m] !NR inverse
    real(kind=WP), parameter :: betas = 0.0008
    real(kind=WP), parameter :: betat = 0.00004
    real(kind=WP)            :: a1, f0, f1, ttmp, qrho

    qrho=betas*qfm-betat*qtm
    ttmp=60.

    if(qrho>0.) then
       ttmp=0.
    else
       do iter=1,5
          a1 = exp(-ttmp*qhw) 
          f0 = 2.* qw * a1 + 9.81 * qrho * ttmp
          f1 =-2.* qw * a1 * qhw + 9.81 * qrho
          ttmp = ttmp - f0 / f1
          ttmp = max(ttmp,10.) 
       enddo
    end if

    obuk=max(ttmp,10.)

end subroutine pmlktmo
!
!
!_______________________________________________________________________________
! calculate non constant background diffusion coefficient as it is also used in 
! KPP mixing scheme (see. oce_ale_mxing_kpp.F90-->subroutine ri_iwmix(viscA, diffK))
! after suggestions from Qiang Wang in FESOM1.4
subroutine Kv0_background_qiang(Kv0_b,lat,dep)
	! Kv0_b ... backround mixing coefficient
	use o_PARAM, only: WP
	implicit none
	real(kind=WP), intent(out) :: Kv0_b
	real(kind=WP), intent(in)  :: lat, dep
	real(kind=WP)              :: aux, ratio
	
	!___________________________________________________________________________
	! set latitudinal and depth dependent background diffusivity after 
	! Q. Wang FESOM1.4 approach
	!_______________________________________________________________________
	aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1.0e-5_WP
	
	!_______________________________________________________________________
	! latitudinal equatorial scaling
	if (abs(lat) < 5.0_WP) then
		ratio = 1.0_WP
	else
		ratio = MIN( 1.0_WP + 9.0_WP * (abs(lat) - 5.0_WP) / 10.0_WP, 10.0_WP )
	end if 
	
	!_______________________________________________________________________
	! latitudinal arctic scaling
	if (lat > 70.0) then
		if (dep <= 50.0)     then
			ratio=4.0_WP + 6.0_WP * (50.0_WP - dep) / 50.0_WP
		else
			ratio=4.0_WP  
		endif   
	end if
	!_______________________________________________________________________
	Kv0_b = aux*ratio
	
end subroutine Kv0_background_qiang
!
!
!_______________________________________________________________________________
! calculate non constant background diffusion coefficient as it is also used in 
! KPP mixing scheme (see. oce_ale_mxing_kpp.F90-->subroutine ri_iwmix(viscA, diffK))
! first implemented my Q.Wang in FESOM1.4
subroutine Kv0_background(Kv0_b,lat,dep)
	! Kv0_b ... backround mixing coefficient
	use o_PARAM, only: WP
	implicit none
	real(kind=WP), intent(out) :: Kv0_b
	real(kind=WP), intent(in)  :: lat, dep
	real(kind=WP)              :: aux, ratio
	
	!___________________________________________________________________________
	! set latitudinal and depth dependent background diffusivity after 
	! Q. Wang FESOM1.4 approach
	!_______________________________________________________________________
	! latitudinal equatorial scaling
	if (abs(lat) < 5.0_WP) then
		ratio = 1.0_WP
	else
		ratio = MIN( 1.0_WP + 9.0_WP * (abs(lat) - 5.0_WP) / 10.0_WP, 10.0_WP )
	end if 
	
	!_______________________________________________________________________
	! latitudinal <70° scaling
	if (lat < 70.0) then
		aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1e-5_WP
	! latitudinal arctic scaling
	else
		aux = (0.6_WP + 1.0598_WP / 3.1415926_WP * ATAN( 4.5e-3_WP * (dep - 2500.0_WP))) * 1.e-6
		ratio=3.0
		if (dep < 80.0)     then
			ratio=1.0
		elseif(dep < 100.0) then
			ratio=2.0  
		endif   
	end if
	!_______________________________________________________________________
	Kv0_b = aux*ratio
	
end subroutine Kv0_background