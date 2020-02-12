MODULE gen_bulk
    ! Compute heat and momentum exchange coefficients
    use o_mesh
    use i_therm_param
    use i_arrays
    use g_forcing_arrays
    use g_forcing_param, only: ncar_bulk_z_wind, ncar_bulk_z_tair, ncar_bulk_z_shum
    use g_parsup
    use o_param, only: WP
    use g_sbf, only: atmdata, i_totfl, i_xwind, i_ywind, i_humi, i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_cloud

    implicit none

    public ncar_ocean_fluxes_mode
    public core_coeff_2z
  
    CONTAINS
    
subroutine ncar_ocean_fluxes_mode 
    ! Compute drag coefficient and the transfer coefficients for evaporation
    ! and sensible heat according to LY2004.
    ! In this routine we assume air temperature and humidity are at the same
    ! height as wind speed. Otherwise, the code should be modified.
    ! There is a parameter z, which sets the height of wind speed. 
    ! For the CORE forcing data, z=10.0 
    !
    ! original note:
    ! Over-ocean fluxes following Large and Yeager (used in NCAR models)           
    ! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
    ! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
    ! Stephen.Griffies@noaa.gov updated the code with the bug fix.  
    ! 
    ! Code from CORE website is adopted to FESOM by Qiang Wang
    ! Reviewed by ??
    !----------------------------------------------------------------------

    integer, parameter :: n_itts = 2
    integer            :: i, j, m
    real(kind=WP) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
    real(kind=WP) :: cd, ce, ch, cd_rt                    ! full drag coefficients @ z
    real(kind=WP) :: zeta_u, zeta_tq, x2, x, psi_m_u, psi_h_u, psi_m_tq, psi_h_tq, stab      ! stability parameters
    real(kind=WP) :: t, ts, q, qs, u, u10, tv, xx, dux, dvy
    real(kind=WP) :: tstar, qstar, ustar, bstar
    real(kind=WP), parameter :: grav = 9.80_WP, vonkarm = 0.40_WP
    real(kind=WP), parameter :: q1=640380._WP, q2=-5107.4_WP    ! for saturated surface specific humidity

    do i=1,myDim_nod2d+eDim_nod2d       
        ! degree celcium to Kelvin
        t      = tair(i) + tmelt 
        ts     = t_oc_array(i) + tmelt  
        
        q      = shum(i)
        qs     = 0.98_WP*q1*inv_rhoair*exp(q2/ts)                   ! L-Y eqn. 5 
        ! virtual potential temperature tv
        tv     = t*(1.0_WP+0.608_WP*q)
        
        ! first guess 10m wind
        dux    = u_wind(i)-u_w(i)
        dvy    = v_wind(i)-v_w(i)
        u      = max(sqrt(dux**2+dvy**2), 0.5_WP)  ! 0.5 m/s floor on wind (undocumented NCAR)
        u10    = u                               
        
        cd_n10 = (2.7_WP/u10+0.142_WP+0.0764_WP*u10)*1.0e-3_WP     ! L-Y eqn. 6a
        cd_n10_rt = sqrt(cd_n10) 
        ce_n10 = 34.6_WP *cd_n10_rt*1.0e-3_WP                      ! L-Y eqn. 6b
        stab   = 0.5_WP + sign(0.5_WP,t-ts)
        ch_n10 = (18.0_WP*stab+32.7_WP*(1.0_WP-stab))*cd_n10_rt*1.e-3_WP ! L-Y eqn. 6c
        
        ! first guess for exchange coeff's at z
        cd     = cd_n10                                 
        ch     = ch_n10
        ce     = ce_n10
        
        !_______________________________________________________________________
        ! Monin-Obukhov iteration
        do j=1,n_itts  
            !___________________________________________________________________
            ! (1) for j==1 calculate initial turbulent scales
            cd_rt = sqrt(cd)
            ustar = cd_rt*u                                        ! L-Y eqn. 7a
            tstar = (ch/cd_rt)*(t-ts)                              ! L-Y eqn. 7b
            qstar = (ce/cd_rt)*(q-qs)                              ! L-Y eqn. 7c
            bstar = grav*(tstar/tv+qstar/(q+1.0_WP/0.608_WP))
            
            !___________________________________________________________________
            ! (2a) calculate stability parameter zeta_u = z_u/L, L...Monin-Obukov length
            zeta_u = vonkarm*bstar*ncar_bulk_z_wind/(ustar*ustar)  ! L-Y eqn. 8a
            zeta_u = sign( min(abs(zeta_u),10.0_WP), zeta_u)       ! undocumented NCAR
            x2     = sqrt(abs(1._WP-16._WP*zeta_u))                ! L-Y eqn. 8b
            x2     = max(x2, 1.0_WP)                               ! undocumented NCAR
            x      = sqrt(x2)
            
            !___________________________________________________________________
            ! calculate integrals of dimensionless flux profiles of momentum 
            ! psi_m... and heat and moisture psi_h...
            if (zeta_u > 0._WP) then
                psi_m_u = -5._WP*zeta_u                            ! L-Y eqn. 8c
                psi_h_u = -5._WP*zeta_u                            ! L-Y eqn. 8c
            else
                psi_m_u = log((1._WP+2._WP*x+x2)*(1.0_WP+x2)/8._WP)-2._WP*(atan(x)-atan(1.0_WP))  ! L-Y eqn. 8d
                psi_h_u = 2._WP*log((1._WP+x2)/2._WP)                                  ! L-Y eqn. 8e
            end if
            
            !___________________________________________________________________
            ! (2b) calculate stability parameter zeta_tq = z_tq/L, L...Monin-Obukov length
            zeta_tq = vonkarm*bstar*ncar_bulk_z_tair/(ustar*ustar)! L-Y eqn. 8a
            zeta_tq = sign( min(abs(zeta_tq),10.0), zeta_tq )     ! undocumented NCAR
            x2      = sqrt(abs(1.-16.*zeta_tq))                   ! L-Y eqn. 8b
            x2      = max(x2, 1.0)                                ! undocumented NCAR
            x       = sqrt(x2)
            
            !___________________________________________________________________
            ! calculate integrals of dimensionless flux profiles of momentum 
            ! psi_m... and heat and moisture psi_h...
            if (zeta_tq > 0.) then
                psi_m_tq = -5.*zeta_tq                             ! L-Y eqn. 8c
                psi_h_tq = -5.*zeta_tq                             ! L-Y eqn. 8c
            else
                psi_m_tq = log((1.+2.*x+x2)*(1+x2)/8.)-2.*(atan(x)-atan(1.0))! L-Y eqn. 8d
                psi_h_tq = 2.*log((1.+x2)/2.)                                ! L-Y eqn. 8e
            end if
            
            !___________________________________________________________________
            ! (3a) shift wind speed to 10m and neutral stability
            u10 = u/(1.0_WP+cd_n10_rt*(log(ncar_bulk_z_wind/10._WP)-psi_m_u)/vonkarm) ! L-Y eqn. 9a !why cd_n10_rt not cd_rt
            
            !___________________________________________________________________
            ! (3b) shift temperature and humidity to wind height
            t = t - tstar/vonkarm*(log(ncar_bulk_z_tair/ncar_bulk_z_wind)+psi_h_u-psi_h_tq)! L-Y eqn. 9b
            q = q - qstar/vonkarm*(log(ncar_bulk_z_shum/ncar_bulk_z_wind)+psi_h_u-psi_h_tq)! L-Y eqn. 9b
            
            !___________________________________________________________________
            ! (4a) update neutral 10m transfer coefficient
            cd_n10 = (2.7_WP/u10+0.142_WP+0.0764_WP*u10)*1.e-3_WP  ! L-Y eqn. 6a again
            cd_n10_rt = sqrt(cd_n10) 
            ce_n10 = 34.6_WP*cd_n10_rt*1.e-3_WP                    ! L-Y eqn. 6b again
            stab = 0.5_WP + sign(0.5_WP,zeta_u)
            ch_n10 = (18.0_WP*stab+32.7_WP*(1.0_WP-stab))*cd_n10_rt*1.e-3_WP ! L-Y eqn. 6c again
            !z0 = 10*exp(-vonkarm/cd_n10_rt)                          ! diagnostic
            
            !_______________________________________________________________________
            ! (4b) shift them to the measurement height z_wind and stability zeta_u
            xx = (log(ncar_bulk_z_wind/10._WP)-psi_m_u)/vonkarm
            cd = cd_n10/(1.0_WP+cd_n10_rt*xx)**2                   ! L-Y 10a
            xx = (log(ncar_bulk_z_wind/10._WP)-psi_h_u)/vonkarm
            ch = ch_n10/(1.0_WP+ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)     ! 10b (corrected code aug2007)
            ce = ce_n10/(1.0_WP+ce_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)     ! 10c (corrected code aug2007)
            
            !_______________________________________________________________________
            ! (5) recompute virtual potential temperature tv and update tubulent 
            !     scales (eqn. 7) at the beginning of next iteration loop
            tv = t*(1.0+0.608*q)
            
        end do
        
        cd_atm_oce_arr(i)=cd
        ch_atm_oce_arr(i)=ch
        ce_atm_oce_arr(i)=ce 
        
    end do

end subroutine ncar_ocean_fluxes_mode
!
!---------------------------------------------------------------------------------------------------
!
subroutine cal_wind_drag_coeff
  ! Compute wind-ice drag coefficient following AOMIP
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------
  
  use o_mesh
  use i_arrays
  use g_forcing_arrays
  use g_parsup
  implicit none

  integer            :: i
  real(kind=WP)      :: ws

  do i=1,myDim_nod2d+eDim_nod2d    
     ws=sqrt(u_wind(i)**2+v_wind(i)**2)
     cd_atm_ice_arr(i)=(1.1_WP+0.04_WP*ws)*1.0e-3_WP
  end do

end subroutine cal_wind_drag_coeff
!
SUBROUTINE nemo_ocean_fluxes_mode
!!----------------------------------------------------------------------
!! ** Purpose : Change model variables according to atm fluxes
!! source of original code: NEMO 3.1.1 + NCAR
!!----------------------------------------------------------------------
   IMPLICIT NONE
   integer             :: i
   real(wp)            :: rtmp    ! temporal real
   real(wp)            :: wndm    ! delta of wind module and ocean curent module
   real(wp)            :: wdx,wdy ! delta of wind x/y and ocean curent x/y
   real(wp)            :: q_sat   ! sea surface specific humidity         [kg/kg]
   real(wp), parameter :: rhoa = 1.22 ! air density
   real(wp), parameter :: cpa  = 1000.5         ! specific heat of air
   real(wp), parameter :: Lv   =    2.5e6       ! latent heat of vaporization
   real(wp), parameter :: Stef =    5.67e-8     ! Stefan Boltzmann constant
   real(wp), parameter :: albo =    0.066       ! ocean albedo assumed to be contant
   real(wp)            :: zst     ! surface temperature in Kelvin
   real(wp)           ::  &
      Cd,       &     ! transfer coefficient for momentum         (tau)
      Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
      Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
      t_zu,     &     ! air temp. shifted at zu                     [K]
      q_zu            ! spec. hum.  shifted at zu               [kg/kg]
   real(wp)           :: zevap, zqsb, zqla, zqlw
!!$OMP PARALLEL
!!$OMP DO
   do i = 1, myDim_nod2D+eDim_nod2d
      wdx  = atmdata(i_xwind,i) - u_w(i) ! wind from data - ocean current ( x direction)
      wdy  = atmdata(i_ywind,i) - v_w(i) ! wind from data - ocean current ( y direction)
      wndm = SQRT( wdx * wdx + wdy * wdy )
      zst  = t_oc_array(i)+273.15_WP

      q_sat = 0.98_WP * 640380._WP / rhoa * EXP( -5107.4_WP / zst )

      call core_coeff_2z(2.0_WP, 10.0_WP, zst, atmdata(i_tair,i), &
                        q_sat, atmdata(i_humi,i), wndm, Cd, Ch, Ce, t_zu, q_zu)
     cd_atm_oce_arr(i)=Cd
     ch_atm_oce_arr(i)=Ch
     ce_atm_oce_arr(i)=Ce
   end do
!!$OMP END DO
!!$OMP END PARALLEL
END SUBROUTINE nemo_ocean_fluxes_mode

!---------------------------------------------------------------------------------------------------
SUBROUTINE core_coeff_2z(zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd, Ch, Ce, T_zu, q_zu)
   !!----------------------------------------------------------------------
   !!                      ***  ROUTINE  core_coeff_2z  ***
   !!
   !! ** Purpose :   Computes turbulent transfert coefficients of surface
   !!                fluxes according to Large & Yeager (2004).
   !!
   !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
   !!      Momentum, Latent and sensible heat exchange coefficients
   !!      Caution: this procedure should only be used in cases when air
   !!      temperature (T_air) and air specific humidity (q_air) are at 2m
   !!      whereas wind (dU) is at 10m.
   !!
   !! References :   Large & Yeager, 2004 : ???
   !! code was adopted from NEMO 3.3.1
   !!----------------------------------------------------------------------
   IMPLICIT NONE
   real(wp)            :: dU10        ! dU                             [m/s]
   real(wp)            :: dT          ! air/sea temperature difference   [K]
   real(wp)            :: dq          ! air/sea humidity difference      [K]
   real(wp)            :: Cd_n10      ! 10m neutral drag coefficient
   real(wp)            :: Ce_n10      ! 10m neutral latent coefficient
   real(wp)            :: Ch_n10      ! 10m neutral sensible coefficient
   real(wp)            :: sqrt_Cd_n10 ! root square of Cd_n10
   real(wp)            :: sqrt_Cd     ! root square of Cd
   real(wp)            :: T_vpot      ! virtual potential temperature    [K]
   real(wp)            :: T_star      ! turbulent scale of tem. fluct.
   real(wp)            :: q_star      ! turbulent humidity of temp. fluct.
   real(wp)            :: U_star      ! turb. scale of velocity fluct.
   real(wp)            :: L           ! Monin-Obukov length              [m]
   real(wp)            :: zeta_u      ! stability parameter at height zu
   real(wp)            :: zeta_t      ! stability parameter at height zt
   real(wp)            :: U_n10       ! neutral wind velocity at 10m     [m]
   real(wp)            :: xlogt , xct , zpsi_hu , zpsi_ht , zpsi_m
   real(wp)            :: stab        ! 1st guess stability test integer
   !!
   real(wp), intent(in)   :: &
      zt,      &     ! height for T_zt and q_zt                   [m]
      zu             ! height for dU                              [m]
   real(wp), intent(in)   ::  &
      sst,      &     ! sea surface temperature              [Kelvin]
      T_zt,     &     ! potential air temperature            [Kelvin]
      q_sat,    &     ! sea surface specific humidity         [kg/kg]
      q_zt,     &     ! specific air humidity                 [kg/kg]
      dU              ! relative wind module |U(zu)-U(0)|       [m/s]
   real(wp), intent(out)  ::  &
      Cd,       &     ! transfer coefficient for momentum         (tau)
      Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
      Ce,       &     ! transfert coefficient for evaporation   (Q_lat)
      T_zu,     &     ! air temp. shifted at zu                     [K]
      q_zu            ! spec. hum.  shifted at zu               [kg/kg]

   integer :: j_itt
   integer,  parameter :: nb_itt = 3   ! number of itterations
   real(wp), parameter ::                        &
   grav   = 9.81_WP,      &  ! gravity
   kappa  = 0.4_WP          ! von Karman's constant
   !!----------------------------------------------------------------------
   !!  * Start
   !! Initial air/sea differences
   dU10 = max(0.5_wp, dU)      !  we don't want to fall under 0.5 m/s
   dT = T_zt - sst
   dq = q_zt - q_sat
   !! Neutral Drag Coefficient :
   stab = 0.5_WP + sign(0.5_wp,dT)                 ! stab = 1  if dT > 0  -> STABLE
   Cd_n10  = 1E-3_WP*( 2.7_WP/dU10 + 0.142_WP + dU10/13.09_WP )
   sqrt_Cd_n10 = sqrt(Cd_n10)
   Ce_n10  = 1E-3_WP*( 34.6_WP * sqrt_Cd_n10 )
   Ch_n10  = 1E-3_WP*sqrt_Cd_n10*(18.0_WP*stab + 32.7_WP*(1.0_WP - stab))
   !! Initializing transf. coeff. with their first guess neutral equivalents :
   Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)
   !! Initializing z_u values with z_t values :
   T_zu = T_zt ;  q_zu = q_zt

   !!  * Now starting iteration loop
   do j_itt=1, nb_itt
      dT = T_zu - sst ;  dq = q_zu - q_sat ! Updating air/sea differences
      T_vpot = T_zu*(1._WP + 0.608_WP*q_zu)      ! Updating virtual potential temperature at zu
      U_star = sqrt_Cd*dU10                ! Updating turbulent scales :   (L & Y eq. (7))
      T_star  = Ch/sqrt_Cd*dT              !
      q_star  = Ce/sqrt_Cd*dq              !
      !!
      L = (U_star*U_star) &                ! Estimate the Monin-Obukov length at height zu
           & / (kappa*grav/T_vpot*(T_star*(1._WP+0.608_WP*q_zu) + 0.608_WP*T_zu*q_star))
      !! Stability parameters :
      zeta_u  = zu/L  ;  zeta_u = sign( min(abs(zeta_u),10.0_WP), zeta_u )
      zeta_t  = zt/L  ;  zeta_t = sign( min(abs(zeta_t),10.0_WP), zeta_t )
      zpsi_hu = psi_h(zeta_u)
      zpsi_ht = psi_h(zeta_t)
      zpsi_m  = psi_m(zeta_u)
      !!
      !! Shifting the wind speed to 10m and neutral stability : (L & Y eq.(9a))
      !   U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u)))
      !   In very rare low-wind conditions, the old way of estimating the
      !   neutral wind speed at 10m leads to a negative value that causes the code
      !   to crash. To prevent this a threshold of 0.25m/s is now imposed.
      U_n10 = max(0.25_WP , dU10/(1._WP + sqrt_Cd_n10/kappa*(log(zu/10._WP) - zpsi_m)))
      !!
      !! Shifting temperature and humidity at zu :          (L & Y eq. (9b-9c))
      !T_zu = T_zt - T_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
      T_zu = T_zt - T_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
      !q_zu = q_zt - q_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
      q_zu = q_zt - q_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
      !!
      !! q_zu cannot have a negative value : forcing 0
      stab = 0.5_WP + sign(0.5_WP,q_zu) ;  q_zu = stab*q_zu
      !!
      !! Updating the neutral 10m transfer coefficients :
      Cd_n10  = 1E-3_WP * (2.7_WP/U_n10 + 0.142_WP + U_n10/13.09_WP)    ! L & Y eq. (6a)
      sqrt_Cd_n10 = sqrt(Cd_n10)
      Ce_n10  = 1E-3_WP * (34.6_WP * sqrt_Cd_n10)                 ! L & Y eq. (6b)
      stab    = 0.5_WP + sign(0.5_wp,zeta_u)
      Ch_n10  = 1E-3_WP*sqrt_Cd_n10*(18._WP*stab + 32.7_WP*(1-stab)) ! L & Y eq. (6c-6d)
      !!
      !!
      !! Shifting the neutral 10m transfer coefficients to (zu,zeta_u) :
      !xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u))
      xct = 1._WP + sqrt_Cd_n10/kappa*(log(zu/10._WP) - zpsi_m)
      Cd = Cd_n10/(xct*xct) ; sqrt_Cd = sqrt(Cd)
      !!
      !xlogt = log(zu/10.) - psi_h(zeta_u)
      xlogt = log(zu/10._WP) - zpsi_hu
      !!
      xct = 1._WP + Ch_n10*xlogt/kappa/sqrt_Cd_n10
      Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
      !!
      xct = 1._WP + Ce_n10*xlogt/kappa/sqrt_Cd_n10
      Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
      !!
         !!
   end do
   !!
END SUBROUTINE core_coeff_2z

FUNCTION psi_h( zta )
   !! Psis, L & Y eq. (8c), (8d), (8e)
   !-------------------------------------------------------------------------------
   real(wp)             :: X2
   real(wp)             :: X
   real(wp)             :: stabit
   !
   real(wp), intent(in) ::   zta
   real(wp)             ::   psi_h
   !-------------------------------------------------------------------------------
   X2 = sqrt(abs(1._WP - 16._WP*zta))  ;  X2 = max(X2 , 1._WP) ;  X  = sqrt(X2)
   stabit    = 0.5_WP + sign(0.5_WP,zta)
   psi_h = -5._WP*zta*stabit  &                                       ! Stable
     &    + (1._WP - stabit)*(2._WP*log( (1._WP + X2)/2._WP ))                 ! Unstable
END FUNCTION psi_h

FUNCTION psi_m( zta )
!! Psis, L & Y eq. (8c), (8d), (8e)
!-------------------------------------------------------------------------------
   real(wp)             :: X2
   real(wp)             :: X
   real(wp)             :: stabit
   !!
   real(wp), intent(in) ::   zta
   real(wp), parameter  :: pi = 3.141592653589793_WP
   real(wp)             :: psi_m
   !-------------------------------------------------------------------------------

   X2 = sqrt(abs(1._WP - 16._WP*zta))  ;  X2 = max(X2 , 1.0_WP) ;  X  = sqrt(X2)
   stabit    = 0.5_WP + sign(0.5_WP,zta)
   psi_m = -5._WP*zta*stabit  &                                                          ! Stable
      &    + (1._WP - stabit)*(2._WP*log((1._WP + X)/2._WP) + log((1._WP + X2)/2._WP) - 2._WP*atan(X) + pi/2._WP)  ! Unstable
   !
END FUNCTION psi_m
END MODULE gen_bulk
