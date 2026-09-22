!> @brief
!> oce_density_kernels.F90
!! Pure scalar equation-of-state kernels: no mesh, no partitioning, no
!! diagnostics -- only WP and S_ref_anomaly from o_PARAM.
!!
!! Keep this module's dependencies at o_PARAM only. g_cvmix_kpp uses it, and
!! oce_ale_pressure_bv depends on `diagnostics`, which depends on g_cvmix_kpp, so
!! a dependency added here on any of those three closes a module cycle.

module oce_density_kernels
    use o_PARAM, only: WP, S_ref_anomaly

    implicit none

    private
    public :: densityJM_components

contains

!===============================================================================
! Computes components of the Jackett-McDougall equation of state for seawater.
! This split-form approach separates density into surface density and pressure
! derivatives, enabling efficient computation of both potential and in-situ 
! density from a single EOS call.
!
! INPUT:
!   t          - in-situ temperature (°C)
!   s          - salinity (psu)
!
! OUTPUT:
!   bulk_0     - density at surface pressure (P=0 dbar), kg/m³
!                This is approximately the potential density ρ₀(T,S,P=0)
!   bulk_pz    - first derivative of density w.r.t. pressure: ∂ρ/∂P, kg/m³/dbar
!   bulk_pz2   - second derivative of density w.r.t. pressure: ∂²ρ/∂P², kg/m³/dbar²
!   rhopot     - potential density referenced to surface (σ₀), kg/m³
!
! USAGE - Computing In-Situ Density:
!   In-situ density at depth z (negative, in meters) is computed via Taylor 
!   expansion around surface pressure:
!
!   ρ_insitu(z) = bulk_0 + z·(bulk_pz + z·bulk_pz2)
!
!   This accounts for compressibility: water gets denser with increasing pressure.
!   A compressibility correction is then applied:
!
!   ρ_insitu(z) = ρ_insitu·rhopot / (ρ_insitu + 0.1·z·state_eq)
!
!   The factor 0.1 is a UNIT CONVERSION from depth [m] to pressure [dbar]:
!     P [dbar] ≈ 0.1 × |z| [m]
!   Example: at z = -2000 m → P ≈ 200 dbar
!   This ensures the pressure-dependent correction term has correct units, since
!   bulk_pz and bulk_pz2 are derivatives w.r.t. pressure in decibars.
!
!   Finally, subtract reference density to get density anomaly:
!
!   density_m_rho0(z) = ρ_insitu(z) - density_ref(z)
!
! WHY SPLIT FORM?
!   - Efficient: One EOS call provides both potential and in-situ density
!   - Accurate: Taylor expansion reduces pressure gradient errors
!   - Flexible: Can compute density at any depth without repeated EOS calls
!
! NOTE: 
!   - Potential density (rhopot) has NO pressure effects → water mass properties
!   - In-situ density includes pressure → used for dynamics and PGF calculations
SUBROUTINE densityJM_components(t, s, bulk_0, bulk_pz, bulk_pz2, rhopot)
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
  ! N. Rakowski 2014 the split form
  !---------------------------------------------------------------------------
  
  real(kind=WP),  intent(IN)            :: t,s
  real(kind=WP),  intent(OUT)           :: bulk_0, bulk_pz, bulk_pz2, rhopot
  real(kind=WP)                         :: s_sqrt, s_abs

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

  s_abs = s + S_ref_anomaly   ! EOS needs absolute salinity (S_ref=0 unless use_salt_anomaly)
  s_sqrt = sqrt(s_abs)

  bulk_0 =  a0      + t*(at   + t*(at2  + t*(at3 + t*at4)))      &
          + s_abs* (as  + t*(ast  + t*(ast2 + t*ast3))               &
               + s_sqrt*(ass  + t*(asst + t*asst2)))

  bulk_pz =  ap  + t*(apt  + t*(apt2 + t*apt3))                  &
                  + s_abs*(aps + t*(apst + t*apst2) + s_sqrt*apss)

  bulk_pz2 = ap2 + t*(ap2t + t*ap2t2)		                 &
                + s_abs *(ap2s + t*(ap2st + t*ap2st2))

  rhopot =  b0 + t*(bt + t*(bt2 + t*(bt3  + t*(bt4  + t*bt5))))	 &
               + s_abs*(bs + t*(bst + t*(bst2 + t*(bst3 + t*bst4)))  &
                  + s_sqrt*(bss + t*(bsst + t*bsst2))            &
                       + s_abs* bss2)
end subroutine densityJM_components

end module oce_density_kernels
