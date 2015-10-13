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

real(kind=WP)         :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s
integer               :: node, nz, elem, elnodes(3) ,i
real(kind=WP)         :: rhopot(nl), bulk_0(nl), bulk_pz(nl), bulk_pz2(nl)
real(kind=WP)         :: bulk_up, bulk_dn

  DO node=1, myDim_nod2D+eDim_nod2D
     DO nz=2,nlevels_nod2d(node)-1
        !if (bvfreq(nz,node) < 0) then
        !   Kv(nz,node)=1.0_WP   
        !else
           dz_inv=1.0_WP/(Z(nz-1)-Z(nz))
           shear = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
	           (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
           shear = shear*dz_inv*dz_inv
	   Kv(nz,node) = shear/(shear+5.*max(bvfreq(nz,node),0.0_8)+1.0e-14)  ! To avoid NaNs at start
        !end if
    END DO	
   END DO

   if (use_ice .and. mo_on) then !stress is only partial!!
      DO node=1, myDim_nod2D+eDim_nod2D
         call mo_length(water_flux(node),heat_flux(node), &         
                        stress_atmoce_x(node),stress_atmoce_y(node), &    
                        u_ice(node),v_ice(node),a_ice(node), &                             
                        dt, mixlength(node))

         DO nz = 2,nlevels_nod2d(node)-1
            mo(nz,node) = 0._WP
            if(abs(zbar(nz)) <= mixlength(node)) mo(nz,node)=modiff    ! Potentialy bad place 
         end DO                                                        ! IF inside the internal cycle
      END DO

      DO elem=1, myDim_elem2D
         elnodes=elem2D_nodes(:,elem)
         DO nz=2,nlevels(elem)-1
            Av(nz,elem)= mix_coeff_PP*sum(Kv(nz,elnodes)**2)/3.0_WP + A_ver  &
                          + sum(mo(nz,elnodes))/3.0 
         END DO
      END DO
      DO node=1,myDim_nod2D+eDim_nod2D
         DO nz=2,nlevels_nod2d(node)-1 
            Kv(nz,node)= mix_coeff_PP*Kv(nz,node)**3 + K_ver &
                 + mo(nz,node)
         END DO
      END DO
      
   else ! .not.  (use_ice .and. mo_on)
      
      DO elem=1, myDim_elem2D
         elnodes=elem2D_nodes(:,elem)
         DO nz=2,nlevels(elem)-1
            Av(nz,elem)= mix_coeff_PP*sum(Kv(nz,elnodes)**2)/3.0_WP+A_ver
         END DO
      END DO
      DO node=1,myDim_nod2D+eDim_nod2D
         DO nz=2,nlevels_nod2d(node)-1 
            Kv(nz,node)= mix_coeff_PP*Kv(nz,node)**3+K_ver
         END DO
      END DO

   end if
end subroutine oce_mixing_pp
! ========================================================================
subroutine mo_length(water_flux,heat_flux,stress_x,stress_y,  &
       u_ice,v_ice,a_ice,dt,mixlength)
    ! vertical mixing scheme of Timmermann and Beckmann, 2004.
    ! computes the mixing length derived from the Monin-Obukhov length
    ! Ralph Timmermann, 14.06.2006

    implicit none

    real*8              :: water_flux, heat_flux, stress_x, stress_y
    real*8              :: u_ice, v_ice, a_ice, uabs
    real*8              :: dt, ret, rtc, mixlength
    real*8              :: qfm, qtm, qw
    real*8              :: ustar,tau, obuk
    real(kind=8), parameter :: cosgam = 0.913632  ! cos(24.*3.14/180.)

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
    implicit none

    integer           :: iter
    real*8            :: qtm, qfm, qw, obuk
    real*8, parameter :: qhw   = 1/7.0                 ! [1/m] !NR inverse
    real*8, parameter :: betas = 0.0008
    real*8, parameter :: betat = 0.00004
    real*8            :: a1, f0, f1, ttmp, qrho

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
