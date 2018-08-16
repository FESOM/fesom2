! Compute heat and momentum exchange coefficients

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
  
  use o_mesh
  use i_therm_param
  use i_arrays
  use g_forcing_arrays
  use g_parsup
  implicit none

  integer, parameter :: n_itts = 2
  integer            :: i, j, m
  real(kind=WP) :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real(kind=WP) :: cd, ce, ch, cd_rt                    ! full drag coefficients @ z
  real(kind=WP) :: zeta, x2, x, psi_m, psi_h, stab      ! stability parameters
  real(kind=WP) :: t, ts, q, qs, u, u10, tv, xx, dux, dvy
  real(kind=WP) :: tstar, qstar, ustar, bstar
  real(kind=WP), parameter :: grav = 9.80, vonkarm = 0.40
  real(kind=WP), parameter :: q1=640380., q2=-5107.4    ! for saturated surface specific humidity
  real(kind=WP), parameter :: zz = 10.0

  do i=1,myDim_nod2d+eDim_nod2d       
     t=tair(i) + tmelt					      ! degree celcium to Kelvin
     ts=t_oc_array(i) + tmelt				      !
     q=shum(i)
     qs=0.98*q1*inv_rhoair*exp(q2/ts) 			      ! L-Y eqn. 5 
     tv = t*(1.0+0.608*q)
     dux=u_wind(i)-u_w(i)
     dvy=v_wind(i)-v_w(i)
     u = max(sqrt(dux**2+dvy**2), 0.5)           	      ! 0.5 m/s floor on wind (undocumented NCAR)
     u10 = u                                                  ! first guess 10m wind

     cd_n10 = (2.7/u10+0.142+0.0764*u10)*1.0e-3                ! L-Y eqn. 6a
     cd_n10_rt = sqrt(cd_n10) 
     ce_n10 = 34.6 *cd_n10_rt*1.0e-3       		      ! L-Y eqn. 6b
     stab = 0.5 + sign(0.5,t-ts)
     ch_n10 = (18.0*stab+32.7*(1.0-stab))*cd_n10_rt*1.e-3      ! L-Y eqn. 6c

     cd = cd_n10                                 	      ! first guess for exchange coeff's at z
     ch = ch_n10
     ce = ce_n10
     do j=1,n_itts                                            ! Monin-Obukhov iteration
        cd_rt = sqrt(cd)
        ustar    = cd_rt*u                                    ! L-Y eqn. 7a
        tstar    = (ch/cd_rt)*(t-ts)              	      ! L-Y eqn. 7b
        qstar    = (ce/cd_rt)*(q-qs)              	      ! L-Y eqn. 7c
        bstar    = grav*(tstar/tv+qstar/(q+1.0/0.608))
        zeta     = vonkarm*bstar*zz/(ustar*ustar) 	      ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta )          ! undocumented NCAR
        x2 = sqrt(abs(1.-16.*zeta))                           ! L-Y eqn. 8b
        x2 = max(x2, 1.0)                                     ! undocumented NCAR
        x = sqrt(x2)

        if (zeta > 0.) then
           psi_m = -5.*zeta                                    ! L-Y eqn. 8c
           psi_h = -5.*zeta                                    ! L-Y eqn. 8c
        else
           psi_m = log((1.+2.*x+x2)*(1+x2)/8.)-2.*(atan(x)-atan(1.0))  ! L-Y eqn. 8d
           psi_h = 2.*log((1.+x2)/2.)                                  ! L-Y eqn. 8e
        end if

        u10 = u/(1.0+cd_n10_rt*(log(zz/10.)-psi_m)/vonkarm)        ! L-Y eqn. 9 !why cd_n10_rt not cd_rt
        cd_n10 = (2.7/u10+0.142+0.0764*u10)*1.e-3                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10) 
        ce_n10 = 34.6*cd_n10_rt*1.e-3                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1.0-stab))*cd_n10_rt*1.e-3       ! L-Y eqn. 6c again
        !z0 = 10*exp(-vonkarm/cd_n10_rt)                          ! diagnostic

        xx = (log(zz/10.)-psi_m)/vonkarm
        cd = cd_n10/(1.0+cd_n10_rt*xx)**2             		  ! L-Y 10a
        xx = (log(zz/10.)-psi_h)/vonkarm
        ch = ch_n10/(1.0+ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)     ! 10b (corrected code aug2007)
        ce = ce_n10/(1.0+ce_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)     ! 10c (corrected code aug2007)
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

  integer            :: i, m
  real(kind=WP)      :: ws

  do i=1,myDim_nod2d+eDim_nod2d    
     ws=sqrt(u_wind(i)**2+v_wind(i)**2)
     cd_atm_ice_arr(i)=(1.1+0.04*ws)*1.0e-3
  end do

end subroutine cal_wind_drag_coeff
!
!---------------------------------------------------------------------------------------------------
