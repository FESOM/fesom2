!=======================================================================
subroutine oce_mixing_pp
    !  Compute Richardson number dependent Av and Kv following
    !  Pacanowski and Philander, 1981
    !  Av = Avmax * factor**2 + Av0,
    !  Kv = Kvmax * factor**3 + Kv0, 
    !  factor=1/(1+10Ri),  Ri=N**2/(dU/dz)**2 is the Richardson number
    !                     N is the buoyancy frequency
    !                     dU/dz is the vertical velocity shear
    !  Avmax, Kvmax are tunable
    
    ! Output: Av(2:nlevels(elem)-1,:)  == vert. visc. coeff. at zbar(2:...)
    !         Kv(2:nlevels_nod2D(node),:) == vert. diff. coeff. at zbar(2:...)   
    
USE o_PARAM
USE o_MESH
USE o_ARRAYS
USE g_PARSUP
USE g_config
use i_arrays
IMPLICIT NONE

real(kind=WP)         :: dz_inv, bv, shear, a, rho_up, rho_dn, t, s
integer               :: node, nz, elem, elnodes(3) ,i
real(kind=WP)         :: wndmix=0.001
real(kind=WP)         :: rhopot(nl), bulk_0(nl), bulk_pz(nl), bulk_pz2(nl)
real(kind=WP)         :: bulk_up, bulk_dn

  DO node=1, myDim_nod2D+eDim_nod2D


     DO nz=1,nlevels_nod2d(node)-1
        t=tr_arr(nz, node,1)
        s=tr_arr(nz, node,2)
        call densityJM_components(t, s, bulk_0(nz), bulk_pz(nz), bulk_pz2(nz), rhopot(nz))
     END DO

     DO nz=2,nlevels_nod2d(node)-1

        bulk_up = bulk_0(nz-1) + zbar(nz)*(bulk_pz(nz-1) + zbar(nz)*bulk_pz2(nz-1)) 
        bulk_dn = bulk_0(nz)   + zbar(nz)*(bulk_pz(nz)   + zbar(nz)*bulk_pz2(nz))

!NR "- density_0" omitted, because "rho_up-rho_dn" is needed, thus, "-density_0"
!NR may only degrade the numerical accuracy of the substraction.

        rho_up = bulk_up*rhopot(nz-1) / (bulk_up + 0.1*zbar(nz)) !- density_0 
        rho_dn = bulk_dn*rhopot(nz)   / (bulk_dn + 0.1*zbar(nz)) !- density_0 

  
        if (rho_up-rho_dn > 0) then
           Kv(nz,node)=1.0_WP   
        else
        ! both densities are referenced to the plane 
        ! where flux is computed (zbar(nz)

!NR compute bv only where needed.
           dz_inv=1.0_WP/(Z(nz-1)-Z(nz))  
           bv  = -g*dz_inv*(rho_up-rho_dn)/density_0
           
           shear = (Unode(1,nz-1,node)-Unode(1,nz,node))**2 +&
	           (Unode(2,nz-1,node)-Unode(2,nz,node))**2 
           shear = shear*dz_inv*dz_inv
	   Kv(nz,node) = shear/(shear+5.*bv+1.0e-14)  ! To avoid NaNs at start
        end if
     END DO	

   END DO

!NR If-statement outside loop. Hopefully not to unreadable.

   if (use_ice .and. mo_on) then !stress is only partial!!
      DO node=1, myDim_nod2D+eDim_nod2D
         !NR mixlength is constant for the column, right?!
         !!! !NR NB! mixlength undefined on first entry, but used!!! 
         call mo_length(water_flux(node),heat_flux(node), &         
                        stress_atmoce_x(node),stress_atmoce_y(node), &    
                        u_ice(node),v_ice(node),a_ice(node), &                             
                        dt, mixlength(node))

         DO nz = 2,nlevels_nod2d(node)-1
            mo(nz,node) = 0._WP
            if(abs(zbar(nz)) <= mixlength(node)) mo(nz,node)=modiff
         end DO
      END DO

      DO elem=1, myDim_elem2D
         elnodes=elem2D_nodes(:,elem)
         DO nz=2,nlevels(elem)-1
            Av(nz,elem)= mix_coeff_PP*sum(Kv(nz,elnodes)**2)/3.0_WP + A_ver  &
                          + sum(mo(nz,elnodes))/3.0 
            !        if (Av(nz,elem)<wndmix) Av(nz,elem)=wndmix 
         END DO
      END DO
      DO node=1,myDim_nod2D+eDim_nod2D
         DO nz=2,nlevels_nod2d(node)-1 
            Kv(nz,node)= mix_coeff_PP*Kv(nz,node)**3 + K_ver &
                 + mo(nz,node)
!        if (Av(nz,node)<wndmix) Av(nz,node)=wndmix 
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


!where (Av(1,:)<wndmix) 
!	Av(:,node)=wndmix
!end where

!where (Kv(1,:)<wndmix) 
!	Kv(:,node)=wndmix
!end where

!Av(1:10, :)=0.1
!Kv(1:10, :)=0.01

!Av=1.e-4
!Kv=1.e-5
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
!===========================================================================
SUBROUTINE densityJM_local(t, s, pz, rho_out)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  ! Ralph Timmermann, August 2005
  !---------------------------------------------------------------------------

  real(kind=WP), intent(IN)  :: t,s,pz
  real(kind=WP), intent(OUT) :: rho_out                 
  real(kind=WP)              :: rhopot, bulk
  real(kind=WP)              :: bulk_0, bulk_pz, bulk_pz2
  !compute secant bulk modulus

  call densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)

  bulk = bulk_0 + pz*(bulk_pz + pz*bulk_pz2) 

  rho_out = bulk*rhopot / (bulk + 0.1*pz) - density_0

end subroutine densityJM_local
		
!===========================================================================
SUBROUTINE densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)
USE o_MESH
USE o_ARRAYS
USE o_PARAM
use g_PARSUP, only: par_ex,pe_status
IMPLICIT NONE

  !
  ! - calculates in-situ density as a function of potential temperature
  !   (relative to the surface)
  !   using the Jackett and McDougall equation of state
  !   (Copyright (c) 1992, CSIRO, Australia)
  ! - has been derived from the SPEM subroutine rhocal
  !
  ! Ralph Timmermann, August 2005
  !---------------------------------------------------------------------------

  real(kind=WP), intent(IN)  :: t,s
  real(kind=WP), intent(OUT) :: bulk_0, bulk_pz, bulk_pz2, rhopot
  real(kind=WP)              :: s_sqrt

  real(kind=WP), parameter   :: a0    = 19092.56,     at   = 209.8925
  real(kind=WP), parameter   :: at2   = -3.041638,    at3  = -1.852732e-3
  real(kind=WP), parameter   :: at4   = -1.361629e-5
  real(kind=WP), parameter   :: as    = 104.4077,     ast  = -6.500517
  real(kind=WP), parameter   :: ast2  = .1553190,     ast3 = 2.326469e-4 
  real(kind=WP), parameter   :: ass   = -5.587545,    asst = 0.7390729 
  real(kind=WP), parameter   :: asst2 = -1.909078e-2
  real(kind=WP), parameter   :: ap    = -4.721788e-1, apt  = -1.028859e-2
  real(kind=WP), parameter   :: apt2  = 2.512549e-4,  apt3 = 5.939910e-7 
  real(kind=WP), parameter   :: aps   = 1.571896e-2,  apst = 2.598241e-4 
  real(kind=WP), parameter   :: apst2 = -7.267926e-6, apss = -2.042967e-3
  real(kind=WP), parameter   :: ap2   = 1.045941e-5,  ap2t = -5.782165e-10 
  real(kind=WP), parameter   :: ap2t2 = 1.296821e-7
  real(kind=WP), parameter   :: ap2s  = -2.595994e-7,ap2st = -1.248266e-9 
  real(kind=WP), parameter   :: ap2st2= -3.508914e-9
  
  real(kind=WP), parameter   :: b0 = 999.842594,    bt  = 6.793952e-2
  real(kind=WP), parameter   :: bt2 = -9.095290e-3, bt3 = 1.001685e-4
  real(kind=WP), parameter   :: bt4 = -1.120083e-6, bt5 = 6.536332e-9
  real(kind=WP), parameter   :: bs = 0.824493,      bst = -4.08990e-3
  real(kind=WP), parameter   :: bst2 = 7.64380e-5,  bst3 = -8.24670e-7		
  real(kind=WP), parameter   :: bst4 = 5.38750e-9	
  real(kind=WP), parameter   :: bss = -5.72466e-3,  bsst = 1.02270e-4
  real(kind=WP), parameter   :: bsst2 = -1.65460e-6,bss2 = 4.8314e-4

  !compute secant bulk modulus


!NR The write statement destroys vectorization.
!NR This checks should be performed in one particular loop, with 
!NR a flag triggering the write-statement outside the loop. 
!NR Good place: directly at or after salinity computation.
        !if (s < 0.) then
        !   write (*,*)'s<0 happens!',t,s
        !   pe_status=1
        !endif

  s_sqrt = sqrt(s)

  bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
          + s* (as  + t*(ast  + t*(ast2 + t*ast3))               &
               + s_sqrt*(ass  + t*(asst + t*asst2)))                

  bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                  + s*(aps + t*(apst + t*apst2) + s_sqrt*apss)

  bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)		                 &
                + s *(ap2s + t*(ap2st + t*ap2st2))

  rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))	 &
               + s*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                  + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                       + s* bss2)

end subroutine densityJM_components
